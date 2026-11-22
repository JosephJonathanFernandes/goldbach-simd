// goldbach.cpp - fast, multithreaded, segmented Goldbach verifier.
//
// Verifies that every even n in [4, N] is the sum of two primes.
//
// Key ideas (same strategy used by large-scale verification projects):
//   1. Segmented sieve: memory is O(sqrt(N) + threads * segment), not O(N).
//      Each segment is an odd-only bitmap (1 bit per odd number) that fits in L2.
//   2. "Small p" search: for every even n the smallest prime p with n-p prime
//      is tiny (< 10^4 for all n verified so far, up to 4*10^18). So we never
//      need the full prime list - only a few hundred small primes plus the
//      sieve window [lo - PMAX, hi] around the evens being checked.
//   3. Bit-parallel check: 64 consecutive evens are tested at once. For a fixed
//      prime p, the primality of (n - p) for n = E, E+2, ..., E+126 is exactly
//      64 consecutive bits of the odd-only bitmap, so one shifted 64-bit load
//      tests 64 numbers. We OR these together until all 64 evens are covered.
//   4. Multithreading: std::thread workers grab chunks from an atomic counter.
//   5. Fallback: if an even number is not covered by small primes, it is
//      re-checked exhaustively with deterministic Miller-Rabin (only a true
//      counterexample would fail this).
//
// Build:  g++ -O3 -march=native -std=c++17 -pthread goldbach.cpp -o goldbach
// Usage:  goldbach [N] [threads]          e.g.  goldbach 1e10 8

#define __USE_MINGW_ANSI_STDIO 1   // correct %llu on old MinGW (no-op elsewhere)
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#if defined(__AVX2__) || defined(__SSE2__)
#  include <immintrin.h>
#endif

typedef uint64_t u64;
typedef uint32_t u32;

// ---------------------------------------------------------------------------
// Minimal threading shim: std::thread where available, otherwise native Win32
// (old mingw.org GCC with the "win32" thread model ships without std::thread).
// ---------------------------------------------------------------------------
#if defined(_WIN32) && defined(__GNUC__) && !defined(_GLIBCXX_HAS_GTHREADS)
#  define GB_WIN32_THREADS 1
#  define WIN32_LEAN_AND_MEAN
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#  undef small
#else
#  include <mutex>
#  include <thread>
#endif

namespace gb {
#ifdef GB_WIN32_THREADS
struct Mutex {
    CRITICAL_SECTION cs;
    Mutex() { InitializeCriticalSection(&cs); }
    ~Mutex() { DeleteCriticalSection(&cs); }
    void lock() { EnterCriticalSection(&cs); }
    void unlock() { LeaveCriticalSection(&cs); }
};
class Thread {
    HANDLE h_ = nullptr;
    struct Ctx { void (*fn)(void*); void* arg; };
    static DWORD WINAPI tramp(LPVOID p) {
        Ctx* c = (Ctx*)p; c->fn(c->arg); delete c; return 0;
    }
public:
    Thread(void (*fn)(void*), void* arg) {
        h_ = CreateThread(nullptr, 0, tramp, new Ctx{fn, arg}, 0, nullptr);
    }
    Thread(Thread&& o) : h_(o.h_) { o.h_ = nullptr; }
    void join() { if (h_) { WaitForSingleObject(h_, INFINITE); CloseHandle(h_); h_ = nullptr; } }
};
inline unsigned hardwareThreads() { SYSTEM_INFO si; GetSystemInfo(&si); return si.dwNumberOfProcessors; }
inline void sleepMs(unsigned ms) { Sleep(ms); }
#else
typedef std::mutex Mutex;
class Thread {
    std::thread t_;
public:
    Thread(void (*fn)(void*), void* arg) : t_(fn, arg) {}
    void join() { t_.join(); }
};
inline unsigned hardwareThreads() { return std::thread::hardware_concurrency(); }
inline void sleepMs(unsigned ms) { std::this_thread::sleep_for(std::chrono::milliseconds(ms)); }
#endif
struct LockGuard {
    Mutex& m;
    explicit LockGuard(Mutex& mm) : m(mm) { m.lock(); }
    ~LockGuard() { m.unlock(); }
};
}  // namespace gb

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
static const u64 PMAX = 1u << 14;           // small primes used for the fast path (16384)
static const u64 SEG_EVENS = 1u << 22;      // evens per work chunk (range of 8M numbers)
static const u64 BATCH_CHUNKS = 8;          // contiguous chunks per task (amortizes 64-bit div)
static const u64 DIRECT_LIMIT = 1u << 20;   // below this, verify with a plain sieve
static const u32 PRESIEVE_MAX = 17;         // primes 3..17 come from a 32 KB periodic pattern
static const u32 PRESIEVE2_MAX = 29;        // primes 19..29 from a second 1.6 KB pattern
static const u32 PRESIEVE3_MAX = 43;        // primes 31..43 from a third 253 KB pattern
static const u64 COV_BLOCK = 2048;          // coverage words per sub-block (16 KB, stays in L1)
static const u32 SIEVE_BLOCK = 1u << 18;    // sieve sub-block in bits (32 KB, L1-sized)
static const u32 BLOCKS_PER_CHUNK = (u32)((SEG_EVENS + PMAX + SIEVE_BLOCK - 1) / SIEVE_BLOCK); // 17
static const size_t MAX_DENSE = 96;         // dense-pass primes must stay < 523 (see record note)

// ---------------------------------------------------------------------------
// Deterministic Miller-Rabin for 64-bit (fallback path only)
// ---------------------------------------------------------------------------
static u64 mulmod(u64 a, u64 b, u64 m) {
#if defined(__SIZEOF_INT128__)
    return (u64)((unsigned __int128)a * b % m);
#else
    u64 r = 0;
    a %= m;
    while (b) {
        if (b & 1) { r += a; if (r >= m || r < a) r -= m; }
        u64 a2 = a + a; if (a2 >= m || a2 < a) a2 -= m; a = a2;
        b >>= 1;
    }
    return r;
#endif
}
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static bool isPrimeMR(u64 n) {
    if (n < 2) return false;
    static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 p : bases) { if (n % p == 0) return n == p; }
    u64 d = n - 1; int s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }
    for (u64 a : bases) {  // these 12 bases are deterministic for all n < 2^64
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = false; break; }
        }
        if (comp) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Simple sieve (base primes up to sqrt(N), small primes, and the direct range)
// ---------------------------------------------------------------------------
static std::vector<u32> simpleSieve(u64 limit, std::vector<uint8_t>* flags = nullptr) {
    std::vector<uint8_t> comp(limit + 1, 0);
    std::vector<u32> primes;
    for (u64 i = 2; i <= limit; ++i) {
        if (!comp[i]) {
            primes.push_back((u32)i);
            for (u64 j = i * i; j <= limit; j += i) comp[j] = 1;
        }
    }
    if (flags) {
        flags->assign(limit + 1, 0);
        for (u32 p : primes) (*flags)[p] = 1;
    }
    return primes;
}

// ---------------------------------------------------------------------------
// Shared state
// ---------------------------------------------------------------------------
struct Shared {
    u64 N;
    u64 start;                       // first even handled by segmented workers
    std::vector<u32> basePrimes;     // odd primes up to sqrt(N)
    std::vector<u32> smallPrimes;    // odd primes up to PMAX
    std::vector<u32> smallOff;       // (PMAX - p - 1)/2: bitmap offset of n - p relative to n's slot
    std::atomic<u64> nextChunk{0};
    u64 numChunks = 0;
    std::atomic<bool> failed{false};
    std::atomic<u64> counterexample{0};
    std::atomic<u64> fallbackCount{0};
    std::atomic<u64> chunksDone{0};
    gb::Mutex recordMutex;
    u64 recordP = 0, recordN = 0;    // largest "minimal p" seen
    std::vector<u64> pattern;        // pre-sieved odd-only bitmap, period patternLen bits
    u32 patternLen = 1;
    std::vector<u64> pattern2;       // second pattern for primes in (PRESIEVE_MAX, PRESIEVE2_MAX]
    u32 pattern2Len = 1;
    std::vector<u64> pattern3;       // third pattern for primes in (PRESIEVE2_MAX, PRESIEVE3_MAX]
    u32 pattern3Len = 1;
    size_t firstSieveIdx = 0;        // index in basePrimes of first prime > PRESIEVE_MAX
    size_t densePrimes = 64;         // primes applied in the vectorized pass A
    std::string ckptFile;            // optional checkpoint file path
    gb::Mutex ckptMutex;
    std::atomic<u64> sieveNs{0}, checkNs{0};
    std::atomic<u64> primeCount{0};  // primes in [DIRECT_LIMIT, N], for the pi(N) self-check
};

static bool loadCheckpoint(Shared& sh) {
    if (sh.ckptFile.empty()) return false;
    FILE* f = std::fopen(sh.ckptFile.c_str(), "r");
    if (!f) return false;
    unsigned long long n = 0, chunk = 0, recP = 0, recN = 0, pi = 0, fb = 0;
    int r = std::fscanf(f, "%llu %llu %llu %llu %llu %llu",
                        &n, &chunk, &recP, &recN, &pi, &fb);
    std::fclose(f);
    if (r == 6 && n == sh.N && chunk <= sh.numChunks) {
        sh.nextChunk.store(chunk);
        sh.chunksDone.store(chunk);
        sh.recordP = recP;
        sh.recordN = recN;
        sh.primeCount.store(pi);
        sh.fallbackCount.store(fb);
        return chunk > 0;
    }
    return false;
}

static void saveCheckpoint(Shared& sh, u64 safeChunk) {
    if (sh.ckptFile.empty()) return;
    gb::LockGuard lk(sh.ckptMutex);
    FILE* f = std::fopen(sh.ckptFile.c_str(), "w");
    if (!f) return;
    std::fprintf(f, "%llu %llu %llu %llu %llu %llu\n",
                 (unsigned long long)sh.N,
                 (unsigned long long)safeChunk,
                 (unsigned long long)sh.recordP,
                 (unsigned long long)sh.recordN,
                 (unsigned long long)sh.primeCount.load(),
                 (unsigned long long)sh.fallbackCount.load());
    std::fclose(f);
}

static inline int popcnt64(u64 x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0; while (x) { x &= x - 1; ++c; } return c;
#endif
}

// Number of set bits with index in [a, b).
static u64 countBits(const u64* w, u64 a, u64 b) {
    if (a >= b) return 0;
    u64 wa = a >> 6, wb = (b - 1) >> 6;
    u64 ma = ~0ULL << (a & 63);
    u64 mb = ~0ULL >> (63 - ((b - 1) & 63));
    if (wa == wb) return popcnt64(w[wa] & ma & mb);
    u64 c = popcnt64(w[wa] & ma) + popcnt64(w[wb] & mb);
    for (u64 i = wa + 1; i < wb; ++i) c += popcnt64(w[i]);
    return c;
}

// Odd-only pattern with multiples of the odd primes in (lo, hi] cleared.
// Bit j <-> odd number 2j+1; period = product of those primes.
static void makePattern(const std::vector<u32>& primes, u32 lo, u32 hi,
                        std::vector<u64>& pat, u32& L) {
    L = 1;
    for (u32 p : primes) if (p > lo && p <= hi) L *= p;
    const u32 totalBits = L + 128;   // padding so extract64(pos < L) never wraps
    pat.assign(totalBits / 64 + 2, ~0ULL);
    for (u32 p : primes) {
        if (p <= lo || p > hi) continue;
        // odd multiples of p: 2j+1 = p*(2t+1)  ->  j = (p-1)/2 + p*t
        for (u32 j = (p - 1) / 2; j < totalBits; j += p)
            pat[j >> 6] &= ~(1ULL << (j & 63));
    }
}

static void buildPattern(Shared& sh) {
    makePattern(sh.basePrimes, 2, PRESIEVE_MAX, sh.pattern, sh.patternLen);
    makePattern(sh.basePrimes, PRESIEVE_MAX, PRESIEVE2_MAX, sh.pattern2, sh.pattern2Len);
    makePattern(sh.basePrimes, PRESIEVE2_MAX, PRESIEVE3_MAX, sh.pattern3, sh.pattern3Len);
    while (sh.firstSieveIdx < sh.basePrimes.size() &&
           sh.basePrimes[sh.firstSieveIdx] <= PRESIEVE3_MAX) ++sh.firstSieveIdx;
}

static inline int ctz64(u64 x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0; while (!(x & 1)) { x >>= 1; ++n; } return n;
#endif
}

// Fast 4-prime fused AVX2 kernel: keeps C[w..w+3] in a YMM register across 4 primes,
// cutting L1 load/store traffic on C[] by 4x.
#if defined(__AVX2__)
static inline __m256i loadShifted256(const u64* src, const __m128i& sr,
                                     const __m128i& sl, unsigned s) {
    __m256i a0 = _mm256_loadu_si256((const __m256i*)src);
    if (s == 0) return a0;
    __m256i a1 = _mm256_loadu_si256((const __m256i*)(src + 1));
    return _mm256_or_si256(_mm256_srl_epi64(a0, sr), _mm256_sll_epi64(a1, sl));
}

static inline void orShifted4Block(u64* __restrict C, const u64* __restrict B,
                                   u64 b0, u64 len, const u32* offs4, bool firstGroup) {
    const u64* s0 = B + (offs4[0] >> 6) + b0; const unsigned sh0 = offs4[0] & 63;
    const u64* s1 = B + (offs4[1] >> 6) + b0; const unsigned sh1 = offs4[1] & 63;
    const u64* s2 = B + (offs4[2] >> 6) + b0; const unsigned sh2 = offs4[2] & 63;
    const u64* s3 = B + (offs4[3] >> 6) + b0; const unsigned sh3 = offs4[3] & 63;
    const __m128i sr0 = _mm_cvtsi32_si128((int)sh0), sl0 = _mm_cvtsi32_si128((int)(64 - sh0));
    const __m128i sr1 = _mm_cvtsi32_si128((int)sh1), sl1 = _mm_cvtsi32_si128((int)(64 - sh1));
    const __m128i sr2 = _mm_cvtsi32_si128((int)sh2), sl2 = _mm_cvtsi32_si128((int)(64 - sh2));
    const __m128i sr3 = _mm_cvtsi32_si128((int)sh3), sl3 = _mm_cvtsi32_si128((int)(64 - sh3));
    u64 w = 0;
    for (; w + 4 <= len; w += 4) {
        __m256i v0 = loadShifted256(s0 + w, sr0, sl0, sh0);
        __m256i v1 = loadShifted256(s1 + w, sr1, sl1, sh1);
        __m256i v2 = loadShifted256(s2 + w, sr2, sl2, sh2);
        __m256i v3 = loadShifted256(s3 + w, sr3, sl3, sh3);
        __m256i v  = _mm256_or_si256(_mm256_or_si256(v0, v1), _mm256_or_si256(v2, v3));
        if (!firstGroup) v = _mm256_or_si256(v, _mm256_loadu_si256((const __m256i*)(C + w)));
        _mm256_storeu_si256((__m256i*)(C + w), v);
    }
    for (; w < len; ++w) {
        u64 v = (sh0 ? (s0[w] >> sh0) | (s0[w + 1] << (64 - sh0)) : s0[w])
              | (sh1 ? (s1[w] >> sh1) | (s1[w + 1] << (64 - sh1)) : s1[w])
              | (sh2 ? (s2[w] >> sh2) | (s2[w + 1] << (64 - sh2)) : s2[w])
              | (sh3 ? (s3[w] >> sh3) | (s3[w + 1] << (64 - sh3)) : s3[w]);
        C[w] = firstGroup ? v : (C[w] | v);
    }
}
#endif

// Fast 256-bit / 128-bit / 64-bit shifted OR of `len` 64-bit words:
//   C[w] |= (src[w] >> s) | (src[w+1] << (64-s))
static inline void orShiftedBlock(u64* __restrict C, const u64* __restrict src,
                                  u64 len, unsigned s) {
    u64 w = 0;
#if defined(__AVX2__)
    if (s == 0) {
        for (; w + 4 <= len; w += 4) {
            __m256i c = _mm256_loadu_si256((const __m256i*)(C + w));
            __m256i a = _mm256_loadu_si256((const __m256i*)(src + w));
            _mm256_storeu_si256((__m256i*)(C + w), _mm256_or_si256(c, a));
        }
    } else {
        const __m128i sr = _mm_cvtsi32_si128((int)s);
        const __m128i sl = _mm_cvtsi32_si128((int)(64 - s));
        for (; w + 4 <= len; w += 4) {
            __m256i c  = _mm256_loadu_si256((const __m256i*)(C + w));
            __m256i a0 = _mm256_loadu_si256((const __m256i*)(src + w));
            __m256i a1 = _mm256_loadu_si256((const __m256i*)(src + w + 1));
            __m256i v  = _mm256_or_si256(_mm256_srl_epi64(a0, sr),
                                         _mm256_sll_epi64(a1, sl));
            _mm256_storeu_si256((__m256i*)(C + w), _mm256_or_si256(c, v));
        }
    }
#elif defined(__SSE2__)
    if (s == 0) {
        for (; w + 2 <= len; w += 2) {
            __m128i c = _mm_loadu_si128((const __m128i*)(C + w));
            __m128i a = _mm_loadu_si128((const __m128i*)(src + w));
            _mm_storeu_si128((__m128i*)(C + w), _mm_or_si128(c, a));
        }
    } else {
        const __m128i sr = _mm_cvtsi32_si128((int)s);
        const __m128i sl = _mm_cvtsi32_si128((int)(64 - s));
        for (; w + 2 <= len; w += 2) {
            __m128i c  = _mm_loadu_si128((const __m128i*)(C + w));
            __m128i a0 = _mm_loadu_si128((const __m128i*)(src + w));
            __m128i a1 = _mm_loadu_si128((const __m128i*)(src + w + 1));
            __m128i v  = _mm_or_si128(_mm_srl_epi64(a0, sr),
                                      _mm_sll_epi64(a1, sl));
            _mm_storeu_si128((__m128i*)(C + w), _mm_or_si128(c, v));
        }
    }
#endif
    if (s == 0) {
        for (; w < len; ++w) C[w] |= src[w];
    } else {
        const unsigned r = 64 - s;
        for (; w < len; ++w) C[w] |= (src[w] >> s) | (src[w + 1] << r);
    }
}

static inline u64 extract64(const u64* w, u32 bit) {
    u32 i = bit >> 6, s = bit & 63;
    return s ? (w[i] >> s) | (w[i + 1] << (64 - s)) : w[i];
}
static inline u64 extract64(const u64* w, u64 bit) {
    u64 i = bit >> 6, s = bit & 63;
    return s ? (w[i] >> s) | (w[i + 1] << (64 - s)) : w[i];
}

// Exhaustive check of a single even n (used when small primes are not enough).
static u64 fullCheck(u64 n) {
    for (u64 p = 3; p <= n / 2; p += 2)
        if (isPrimeMR(p) && isPrimeMR(n - p)) return p;
    return 0;
}

static void worker(void* arg) {
    Shared* sh = static_cast<Shared*>(arg);
    // Odd-only bitmap over window [W, hi): bit i <-> odd number W + 2i + 1.
    const u64 windowBits = (SEG_EVENS + PMAX) + 128;
    std::vector<u64> bits(windowBits / 64 + 4);
    u64 localRecordP = 0, localRecordN = 0;
    std::vector<u64> cov(COV_BLOCK);

    // Split basePrimes into:
    //   1. medium (p < SIEVE_BLOCK = 262144): hit every 32 KB L1 block >= 1 time.
    //   2. large  (p >= SIEVE_BLOCK): bucket-sieved per 32 KB L1 sub-block so even
    //      large-prime hits are applied while the 32 KB sub-block is hot in L1 cache.
    size_t medEnd = sh->firstSieveIdx;
    while (medEnd < sh->basePrimes.size() && sh->basePrimes[medEnd] < SIEVE_BLOCK) ++medEnd;
    const size_t nMed = medEnd - sh->firstSieveIdx;

    struct LargeHit { u32 localBit; };
    std::vector<u32> nextIdx(nMed + 1), primeStep(nMed + 1), overlapIdx(nMed + 1);
    std::vector<LargeHit> buckets[BATCH_CHUNKS];
    u64 sieveNs = 0, checkNs = 0, localPrimes = 0;

    for (;;) {
        if (sh->failed.load(std::memory_order_relaxed)) return;
        const u64 cStart = sh->nextChunk.fetch_add(BATCH_CHUNKS);
        if (cStart >= sh->numChunks) break;
        const u64 cEnd = std::min(cStart + BATCH_CHUNKS, sh->numChunks);
        const u64 nBatchChunks = cEnd - cStart;

        // Initialize prime offsets ONCE at the start of this contiguous batch of chunks
        // (eliminates 64-bit divisions on the remaining chunks in the batch).
        const u64 lo0 = sh->start + cStart * 2 * SEG_EVENS;
        const u64 W0 = lo0 - PMAX;
        const u64 hiBatch = std::min(sh->start + cEnd * 2 * SEG_EVENS, sh->N + 2);
        const u64 totalBatchBits = (hiBatch - W0) / 2;

        size_t nMedActive = 0;
        for (size_t k = sh->firstSieveIdx; k < medEnd; ++k) {
            const u64 p = sh->basePrimes[k];
            if (p * p >= hiBatch) break;
            u64 m = std::max(p * p, ((W0 + p) / p) * p);
            if ((m & 1) == 0) m += p;
            nextIdx[nMedActive] = (u32)((m - W0 - 1) / 2);
            primeStep[nMedActive] = (u32)p;
            ++nMedActive;
        }
        // Populate per-chunk buckets for large primes (p >= SIEVE_BLOCK).
        // Each hit is pushed directly into the chunk(s) it lands in, including the
        // PMAX/2 overlap region [0, PMAX/2) shared by chunk b and chunk b-1.
        for (u64 b = 0; b < nBatchChunks; ++b) buckets[b].clear();
        const u64 overlapBits = PMAX / 2;
        for (size_t k = medEnd; k < sh->basePrimes.size(); ++k) {
            const u64 p = sh->basePrimes[k];
            if (p * p >= hiBatch) break;
            u64 m = std::max(p * p, ((W0 + p) / p) * p);
            if ((m & 1) == 0) m += p;
            for (u64 idx = (m - W0 - 1) / 2; idx < totalBatchBits; idx += p) {
                u64 b = idx >> 22;           // SEG_EVENS == 1 << 22
                u64 rem = idx & (SEG_EVENS - 1);
                if (b < nBatchChunks)
                    buckets[b].push_back(LargeHit{(u32)rem});
                if (rem < overlapBits && b > 0 && b - 1 < nBatchChunks)
                    buckets[b - 1].push_back(LargeHit{(u32)(rem + SEG_EVENS)});
            }
        }

        for (u64 c = cStart; c < cEnd; ++c) {
            if (sh->failed.load(std::memory_order_relaxed)) return;

            const u64 lo = sh->start + c * 2 * SEG_EVENS;               // first even
            const u64 hi = std::min(lo + 2 * SEG_EVENS, sh->N + 2);     // exclusive
            const u64 W = lo - PMAX;                                     // even
            const u64 top = hi;                                          // sieve odds < top
            const u64 nbits = (top - W) / 2;
            const u64 nwords = (nbits + 63) / 64 + 2;

            // --- 1. Sieve the window [W, top) (odd numbers only) ---------------
            // Cache-blocked: the window is processed in SIEVE_BLOCK-bit (32 KB) pieces
            // so every crossing pass hits L1; each prime's next index is carried over.
            auto tS = std::chrono::steady_clock::now();
            u64* Bw = bits.data();
            const u32 nb32 = (u32)nbits;
            const u32 stride = (u32)SEG_EVENS;
            const u64* P = sh->pattern.data();
            const u32 L = sh->patternLen;
            const u64* P2 = sh->pattern2.data();
            const u32 L2 = sh->pattern2Len;
            const u64* P3 = sh->pattern3.data();
            const u32 L3 = sh->pattern3Len;
            u32 pos = (u32)((W / 2) % L);
            u32 pos2 = (u32)((W / 2) % L2);
            u32 pos3 = (u32)((W / 2) % L3);
            for (u32 blk = 0; blk < nwords * 64; blk += SIEVE_BLOCK) {
                // 1a. Pre-sieve: copy the periodic patterns for primes 3..PRESIEVE3_MAX.
                //     Window bit i <-> absolute odd index j = W/2 + i.
                const u32 w0 = blk / 64;
                const u32 w1 = (u32)std::min<u64>(nwords, (blk + SIEVE_BLOCK) / 64);
                for (u32 w = w0; w < w1; ++w) {
                    Bw[w] = extract64(P, pos) & extract64(P2, pos2) & extract64(P3, pos3);
                    pos += 64; if (pos >= L) pos -= L;
                    pos2 += 64; if (pos2 >= L2) pos2 -= L2;
                    pos3 += 64; if (pos3 >= L3) pos3 -= L3;
                }
                // 1b. Cross off medium base primes inside this block.
                //     Save nextIdx for the overlap snap if this block crosses `stride`.
                static const uint8_t CLR[8] = {0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};
                uint8_t* Bb = reinterpret_cast<uint8_t*>(Bw);
                const u32 end = std::min<u32>(nb32, blk + SIEVE_BLOCK);
                if (end >= stride && blk < stride) {
                    // Snapshot the first index >= stride for the next chunk's overlap window [0, PMAX/2)
                    for (size_t k = 0; k < nMedActive; ++k) {
                        u32 i = nextIdx[k];
                        const u32 step = primeStep[k];
                        for (; i < stride; i += step) Bb[i >> 3] &= CLR[i & 7];
                        overlapIdx[k] = i - stride;
                        for (; i < end; i += step) Bb[i >> 3] &= CLR[i & 7];
                        nextIdx[k] = i;
                    }
                } else {
                    for (size_t k = 0; k < nMedActive; ++k) {
                        u32 i = nextIdx[k];
                        const u32 step = primeStep[k];
                        const u32 step8 = step << 3;
                        if (i + step8 <= end) {
                            // 8x constant-mask wheel: (i + 8*step) has the same bit-in-byte (r)
                            // and advances the byte index by exactly `step` — zero shifts in loop.
                            u32 b = i >> 3;
                            const u32 r = i & 7;
                            const u32 r1 = r + step, r2 = r1 + step, r3 = r2 + step, r4 = r3 + step,
                                      r5 = r4 + step, r6 = r5 + step, r7 = r6 + step;
                            const u32 d1 = r1 >> 3, d2 = r2 >> 3, d3 = r3 >> 3, d4 = r4 >> 3,
                                      d5 = r5 >> 3, d6 = r6 >> 3, d7 = r7 >> 3;
                            const uint8_t m0 = CLR[r],      m1 = CLR[r1 & 7], m2 = CLR[r2 & 7], m3 = CLR[r3 & 7],
                                          m4 = CLR[r4 & 7], m5 = CLR[r5 & 7], m6 = CLR[r6 & 7], m7 = CLR[r7 & 7];
                            const u32 bEnd = (end - r - step8) >> 3;
                            for (; b <= bEnd; b += step) {
                                Bb[b]      &= m0;
                                Bb[b + d1] &= m1;
                                Bb[b + d2] &= m2;
                                Bb[b + d3] &= m3;
                                Bb[b + d4] &= m4;
                                Bb[b + d5] &= m5;
                                Bb[b + d6] &= m6;
                                Bb[b + d7] &= m7;
                            }
                            i = (b << 3) | r;
                        }
                        for (; i < end; i += step) Bb[i >> 3] &= CLR[i & 7];
                        nextIdx[k] = i;
                    }
                }
            }
            if (c + 1 < cEnd) {
                for (size_t k = 0; k < nMedActive; ++k) nextIdx[k] = overlapIdx[k];
            }
            // 1c. Apply bucketed large primes (p >= SIEVE_BLOCK) that actually hit this chunk.
            {
                static const uint8_t CLR[8] = {0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};
                uint8_t* Bb = reinterpret_cast<uint8_t*>(Bw);
                for (const LargeHit& h : buckets[c - cStart]) {
                    if (h.localBit < nb32)
                        Bb[h.localBit >> 3] &= CLR[h.localBit & 7];
                }
            }
            auto tC = std::chrono::steady_clock::now();
            sieveNs += (u64)std::chrono::duration_cast<std::chrono::nanoseconds>(tC - tS).count();
            // Count primes p in (lo, min(hi, N)) : bits [(lo-W)/2, (min(hi,N)-W)/2)
            localPrimes += countBits(Bw, (lo - W) / 2, (std::min(hi, sh->N) - W) / 2);

            // --- 2. Bit-parallel Goldbach check -----------------------------------
            // Word w of 'cov' covers the 64 evens n = lo + 128w + 2k (k = 0..63);
            // bit k set  <=>  some prime p tried so far has n - p prime.
            // For prime p, (n - p) maps to bitmap bit off_p + 64w + k with
            // off_p = (lo - p - W - 1)/2 = (PMAX - p - 1)/2, independent of w. So the
            // shift is constant per prime and the loop over w is branch-free and SIMD-able.
            const u64* B = bits.data();
            const u64 nEven = (hi - lo) / 2;
            const u64 nCov = (nEven + 63) / 64;
            const size_t K = std::min(sh->densePrimes, sh->smallPrimes.size());
            const u32* offs = sh->smallOff.data();
            for (u64 b0 = 0; b0 < nCov; b0 += COV_BLOCK) {
                const u64 b1 = std::min(nCov, b0 + COV_BLOCK);
                u64* C = cov.data();
                const u64 len = b1 - b0;
                // Pass A: first K primes over the whole sub-block (4-prime fused AVX2 / fallback).
                size_t j = 0;
#if defined(__AVX2__)
                for (; j + 4 <= K; j += 4)
                    orShifted4Block(C, B, b0, len, offs + j, j == 0);
#endif
                if (j == 0) std::fill(C, C + len, 0ULL);
                for (; j < K; ++j) {
                    const u32 off = offs[j];
                    orShiftedBlock(C, B + (off >> 6) + b0, len, off & 63);
                }
                // Pass B: finish the (few) words not yet fully covered, one prime at a time.
                for (u64 w = 0; w < len; ++w) {
                    const u64 E = lo + 128 * (b0 + w);
                    const u64 count = std::min<u64>(64, (hi - E) / 2);
                    const u64 want = (count == 64) ? ~0ULL : ((1ULL << count) - 1);
                    u64 found = C[w];
                    if ((found & want) == want) continue;
                    u64 lastP = 0;
                    const u32 wbase = (u32)(64 * (b0 + w));
                    const u32* offs = sh->smallOff.data();
                    const size_t nsp = sh->smallPrimes.size();
                    for (size_t j = K; j < nsp; ++j) {
                        found |= extract64(B, wbase + offs[j]);
                        if ((found & want) == want) { lastP = sh->smallPrimes[j]; break; }
                    }
                    if ((found & want) != want) {
                        // Rare/never: some even here needs p > PMAX. Check exhaustively.
                        u64 miss = ~found & want;
                        while (miss) {
                            int k = ctz64(miss);
                            miss &= miss - 1;
                            u64 n = E + 2 * (u64)k;
                            sh->fallbackCount.fetch_add(1);
                            u64 p = fullCheck(n);
                            if (!p) {
                                sh->failed = true;
                                sh->counterexample = n;
                                return;
                            }
                            if (p > localRecordP) { localRecordP = p; localRecordN = n; }
                        }
                    } else if (lastP > localRecordP) {
                        // lastP is the largest minimal-p among these 64 evens; find which n.
                        for (u64 k = 0; k < count; ++k) {
                            u64 n = E + 2 * k, best = 0;
                            for (u32 q : sh->smallPrimes) {
                                u64 bb = (n - q - W - 1) / 2;
                                if ((B[bb >> 6] >> (bb & 63)) & 1) { best = q; break; }
                            }
                            if (best == lastP) { localRecordP = lastP; localRecordN = n; break; }
                        }
                    }
                    // Words resolved in pass A have min-p <= smallPrimes[K-1] (< 523, the
                    // record already set in the direct range), so they cannot be records.
                }
            }
            checkNs += (u64)std::chrono::duration_cast<std::chrono::nanoseconds>(
                           std::chrono::steady_clock::now() - tC).count();
            sh->chunksDone.fetch_add(1);
        }
        // Flush thread-local stats after each batch so checkpoints stay accurate.
        sh->sieveNs.fetch_add(sieveNs); sieveNs = 0;
        sh->checkNs.fetch_add(checkNs); checkNs = 0;
        sh->primeCount.fetch_add(localPrimes); localPrimes = 0;
        {
            gb::LockGuard lk(sh->recordMutex);
            if (localRecordP > sh->recordP) { sh->recordP = localRecordP; sh->recordN = localRecordN; }
        }
    }
}

static u64 parseNum(const char* s) {
    // Accept "1000000", "1e9", "4e10" etc.
    double d = std::strtod(s, nullptr);
    return (u64)std::llround(d);
}

int main(int argc, char** argv) {
    using Clock = std::chrono::steady_clock;
    u64 N = argc > 1 ? parseNum(argv[1]) : 1000000000ULL;
    unsigned T = argc > 2 ? (unsigned)std::atoi(argv[2]) : gb::hardwareThreads();
    if (T == 0) T = 1;
    if (N & 1) --N;
    if (N < 4) { std::printf("Nothing to verify.\n"); return 0; }

    auto t0 = Clock::now();
    Shared sh;
    sh.N = N;
    if (argc > 3 && std::isdigit((unsigned char)argv[3][0]))
        sh.densePrimes = std::min<size_t>(MAX_DENSE, (size_t)std::atoi(argv[3]));
    else if (argc > 3)
        sh.ckptFile = argv[3];
    if (argc > 4)
        sh.ckptFile = argv[4];

    // --- Direct verification of small evens [4, min(N, DIRECT_LIMIT)) --------
    const u64 directEnd = std::min(N + 2, DIRECT_LIMIT);  // exclusive
    std::vector<uint8_t> isP;
    std::vector<u32> directPrimes = simpleSieve(directEnd, &isP);
    for (u64 n = 4; n < directEnd; n += 2) {
        bool ok = false;
        for (u32 p : directPrimes) {
            if (p > n / 2) break;
            if (isP[n - p]) {
                ok = true;
                if (p > sh.recordP) { sh.recordP = p; sh.recordN = n; }
                break;
            }
        }
        if (!ok) { std::printf("COUNTEREXAMPLE FOUND: %llu\n", (unsigned long long)n); return 1; }
    }

    // --- Segmented, multithreaded verification of [DIRECT_LIMIT, N] -----------
    if (N >= DIRECT_LIMIT) {
        u64 sq = (u64)std::sqrt((double)N) + 2;
        while (sq * sq > N + 2) --sq;
        for (u32 p : simpleSieve(sq + 1)) if (p > 2) sh.basePrimes.push_back(p);
        for (u32 p : simpleSieve(PMAX)) if (p > 2) sh.smallPrimes.push_back(p);
        for (u32 p : sh.smallPrimes) sh.smallOff.push_back((u32)((PMAX - p - 1) / 2));
        buildPattern(sh);
        sh.start = DIRECT_LIMIT;
        sh.numChunks = ((N + 2 - sh.start) / 2 + SEG_EVENS - 1) / SEG_EVENS;
        if (loadCheckpoint(sh)) {
            std::printf("Resumed from checkpoint '%s' at chunk %llu / %llu (%.1f%%)\n",
                        sh.ckptFile.c_str(),
                        (unsigned long long)sh.chunksDone.load(),
                        (unsigned long long)sh.numChunks,
                        100.0 * sh.chunksDone.load() / sh.numChunks);
        }

        std::vector<gb::Thread> pool;
        pool.reserve(T);
        for (unsigned i = 0; i < T; ++i) pool.emplace_back(worker, (void*)&sh);

        // Progress reporting (poll often, print twice a second)
        unsigned ticks = 0;
        while (sh.chunksDone.load() < sh.numChunks && !sh.failed.load()) {
            gb::sleepMs(10);
            if (++ticks % 50) continue;
            double el = std::chrono::duration<double>(Clock::now() - t0).count();
            u64 done = sh.chunksDone.load();
            double frac = (double)done / sh.numChunks;
            std::fprintf(stderr, "\r  %5.1f%%  %.1fs elapsed", 100.0 * frac, el);
            if (ticks % 200 == 0 && done >= T * BATCH_CHUNKS) {
                // Conservative lower bound on contiguous chunks finished by all threads
                u64 safe = ((done - T * BATCH_CHUNKS) / BATCH_CHUNKS) * BATCH_CHUNKS;
                saveCheckpoint(sh, safe);
            }
        }
        for (auto& th : pool) th.join();
        if (!sh.failed.load()) saveCheckpoint(sh, sh.numChunks);
        std::fprintf(stderr, "\r                                   \r");
    }

    double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    if (sh.failed) {
        std::printf("COUNTEREXAMPLE FOUND: %llu\n", (unsigned long long)sh.counterexample.load());
        return 1;
    }
    // pi(N) cross-check: primes below DIRECT_LIMIT (incl. 2) + primes counted by workers.
    u64 piN = sh.primeCount.load();
    for (u32 p : directPrimes) if (p <= N) ++piN;
    std::printf("Goldbach verified for all even n in [4, %llu]\n", (unsigned long long)N);
    std::printf("  threads           : %u\n", T);
    std::printf("  time              : %.3f s\n", secs);
    std::printf("  throughput        : %.1f M evens/s\n", (N / 2.0) / secs / 1e6);
    std::printf("  largest min p     : %llu  (n = %llu = %llu + %llu)\n",
                (unsigned long long)sh.recordP, (unsigned long long)sh.recordN,
                (unsigned long long)sh.recordP,
                (unsigned long long)(sh.recordN - sh.recordP));
    std::printf("  fallback checks   : %llu\n", (unsigned long long)sh.fallbackCount.load());
    double s = sh.sieveNs.load() / 1e9, c = sh.checkNs.load() / 1e9;
    if (s + c > 0)
        std::printf("  CPU split         : sieve %.2f s (%.0f%%), check %.2f s (%.0f%%)\n",
                    s, 100 * s / (s + c), c, 100 * c / (s + c));
    std::printf("  pi(N)             : %llu", (unsigned long long)piN);
    static const struct { u64 n, pi; } known[] = {
        {1000000ULL, 78498ULL}, {10000000ULL, 664579ULL}, {100000000ULL, 5761455ULL},
        {1000000000ULL, 50847534ULL}, {10000000000ULL, 455052511ULL},
        {100000000000ULL, 4118054813ULL}, {1000000000000ULL, 37607912018ULL}};
    for (auto& k : known)
        if (k.n == N) {
            std::printf("  [%s vs known %llu]", piN == k.pi ? "OK" : "MISMATCH!",
                        (unsigned long long)k.pi);
            if (piN != k.pi) { std::printf("\n"); return 2; }
        }
    std::printf("\n");
    return 0;
}
