// research_suite.cpp - High-throughput AVX2 Computational Number Theory Research Suite
//
// Four original / frontier experiments built on the goldbach-simd AVX2 segmented sieve:
//   1. lemoine : Verify Lemoine's Conjecture (odd n = p + 2q) past the 10^10 literature
//                benchmark up to 10^12+, recording maximal minimal-q records.
//   2. twin    : Twin-Prime Goldbach (even n = p + q, BOTH p and q twin primes) up to
//                10^11+, mapping all exceptions and maximal minimal-twin-prime records.
//   3. carry   : Goldbach Binary Carry Excess (hardware/CS meets number theory):
//                carries(p, q) = popcount(p) + popcount(q) - popcount(p + q) >= v2(n).
//                Measures Delta_c(n) = min_p carries(p, n-p) - v2(n), zero-excess rate,
//                and finds all exceptions where Delta_c(n) > 0.
//   4. entropy : Information-theoretic analysis of the minimal Goldbach witness stream
//                p_min(n) across billions of evens: Shannon entropy H(p_min), lag-1
//                mutual information I(p_min(n); p_min(n+2)), longest constant runs,
//                and mod-6 residue structure.
//
// Build: g++ -O3 -march=native -std=c++14 -Wall -Wno-format research_suite.cpp -o research_suite
// Usage: research_suite <mode: lemoine|twin|carry|entropy> [N] [threads]

#define __USE_MINGW_ANSI_STDIO 1
#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
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
// Portable threading shim (std::thread or Win32 on old MinGW)
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
    static DWORD WINAPI tramp(LPVOID p) { Ctx* c = (Ctx*)p; c->fn(c->arg); delete c; return 0; }
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
// Bit & primality utilities
// ---------------------------------------------------------------------------
static inline int popcnt64(u64 x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0; while (x) { x &= x - 1; ++c; } return c;
#endif
}
static inline int ctz64(u64 x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0; while (!(x & 1)) { x >>= 1; ++n; } return n;
#endif
}
static u64 mulmod(u64 a, u64 b, u64 m) {
#if defined(__SIZEOF_INT128__)
    return (u64)((unsigned __int128)a * b % m);
#else
    u64 r = 0; a %= m;
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
    for (u64 a : bases) {
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

static inline u64 extract64(const u64* w, u32 bit) {
    u32 i = bit >> 6, s = bit & 63;
    return s ? (w[i] >> s) | (w[i + 1] << (64 - s)) : w[i];
}

// ---------------------------------------------------------------------------
// AVX2 4-prime fused shifted OR kernel
// ---------------------------------------------------------------------------
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
            __m256i v  = _mm256_or_si256(_mm256_srl_epi64(a0, sr), _mm256_sll_epi64(a1, sl));
            _mm256_storeu_si256((__m256i*)(C + w), _mm256_or_si256(c, v));
        }
    }
#endif
    if (s == 0) { for (; w < len; ++w) C[w] |= src[w]; }
    else { const unsigned r = 64 - s; for (; w < len; ++w) C[w] |= (src[w] >> s) | (src[w + 1] << r); }
}

// ---------------------------------------------------------------------------
// Shared Sieve Engine (AVX2 5-pattern L1 pre-sieve + 8x wheel + bucket sieve)
// ---------------------------------------------------------------------------
static const u64 SEG_EVENS    = 1u << 22;   // 4M odd numbers per chunk (8M span)
static const u64 BATCH_CHUNKS = 8;
static const u64 DIRECT_LIMIT = 1u << 20;   // 1,048,576
static const u32 PRESIEVE_MAX = 47;
static const u32 SIEVE_BLOCK  = 1u << 18;   // 32 KB L1 sub-block (262,144 bits)
static const u64 COV_BLOCK    = 2048;

enum Mode { MODE_LEMOINE, MODE_TWIN, MODE_CARRY, MODE_ENTROPY };

struct SieveEnv {
    Mode mode;
    u64 N;
    u64 PMAX;           // overlap window in integers (must be multiple of 16384)
    u64 start;
    u64 numChunks = 0;
    std::vector<u32> basePrimes;
    std::vector<u32> smallPrimes; // witnesses (primes for lemoine/entropy, twin primes for twin)
    std::vector<u32> smallOff;
    std::vector<uint8_t> pat[5];
    u32 patLen[5] = {1, 1, 1, 1, 1};
    size_t firstSieveIdx = 0;
    size_t densePrimes = 64;

    std::atomic<u64> nextChunk{0};
    std::atomic<u64> chunksDone{0};
    std::atomic<bool> failed{false};
    std::atomic<u64> counterexample{0};
    std::atomic<u64> fallbackCount{0};

    gb::Mutex statMutex;
    u64 recordP = 0, recordN = 0;
    std::vector<u64> exceptions; // for twin-prime Goldbach or carry excess > 0

    // Carry-excess stats
    u64 carryTotal = 0;
    u64 carryHist[64] = {0};          // histogram of Delta_c(n)
    u64 carryByPop[64] = {0};         // count of Delta_c(n) > 0 by popcount(n)
    u64 totalByPop[64] = {0};
    u64 maxCarryExcess = 0, maxCarryN = 0;

    // Entropy / information-theoretic stats on p_min(n)
    static const size_t HIST_BINS = 256; // first 255 odd primes + overflow bin
    u64 histAll[HIST_BINS] = {0};
    u64 histMod6[3][HIST_BINS] = {{0}};
    static const size_t TRANS_K = 32;    // lag-1 transition matrix for top 31 primes + other
    u64 trans[TRANS_K][TRANS_K] = {{0}};
    u64 maxRunLen = 0, maxRunPrime = 0, maxRunStart = 0;
    u64 totalEvens = 0;
};

static void makeBytePattern(const std::vector<u32>& primes, u32 lo, u32 hi,
                            std::vector<uint8_t>& pat, u32& L) {
    L = 1;
    for (u32 p : primes) if (p > lo && p <= hi) L *= p;
    const u32 totalBytes = L + 64;
    const u32 totalBits = totalBytes * 8;
    pat.assign(totalBytes, 0xFF);
    for (u32 p : primes) {
        if (p <= lo || p > hi) continue;
        for (u32 j = (p - 1) / 2; j < totalBits; j += p)
            pat[j >> 3] &= (uint8_t)~(1u << (j & 7));
    }
}

static void buildPatterns(SieveEnv& sh) {
    static const u32 bounds[6] = {2, 13, 23, 31, 41, PRESIEVE_MAX};
    for (int i = 0; i < 5; ++i)
        makeBytePattern(sh.basePrimes, bounds[i], bounds[i + 1], sh.pat[i], sh.patLen[i]);
    while (sh.firstSieveIdx < sh.basePrimes.size() &&
           sh.basePrimes[sh.firstSieveIdx] <= PRESIEVE_MAX) ++sh.firstSieveIdx;
}

// Extract twin-prime bitmap T from prime bitmap B:
// T[i] = 1 iff B[i] == 1 and (B[i-1] == 1 or B[i+1] == 1).
static void buildTwinBitmap(const u64* B, u64* T, u64 nwords) {
    u64 prev = 0;
    for (u64 w = 0; w < nwords; ++w) {
        u64 cur = B[w];
        u64 next = (w + 1 < nwords) ? B[w + 1] : 0;
        u64 shL = (cur << 1) | (prev >> 63); // bit i gets B[i-1]
        u64 shR = (cur >> 1) | (next << 63); // bit i gets B[i+1]
        T[w] = cur & (shL | shR);
        prev = cur;
    }
}

static bool isTwinPrimeMR(u64 x) {
    return isPrimeMR(x) && (isPrimeMR(x - 2) || isPrimeMR(x + 2));
}

static void worker(void* arg) {
    SieveEnv* sh = static_cast<SieveEnv*>(arg);
    const u64 PMAX = sh->PMAX;
    const u64 windowBits = (SEG_EVENS + PMAX / 2) + 256;
    std::vector<u64> bits(windowBits / 64 + 8);
    std::vector<u64> twinBits;
    if (sh->mode == MODE_TWIN) twinBits.resize(bits.size());
    std::vector<u64> cov(COV_BLOCK);

    size_t medEnd = sh->firstSieveIdx;
    while (medEnd < sh->basePrimes.size() && sh->basePrimes[medEnd] < SIEVE_BLOCK) ++medEnd;
    const size_t nMed = medEnd - sh->firstSieveIdx;

    struct LargeHit { u32 localBit; };
    std::vector<u32> nextIdx(nMed + 1), primeStep(nMed + 1), overlapIdx(nMed + 1);
    std::vector<LargeHit> buckets[BATCH_CHUNKS];

    u64 localRecordP = 0, localRecordN = 0;
    std::vector<u64> localExceptions;
    u64 lCarryTotal = 0, lCarryHist[64] = {0}, lCarryByPop[64] = {0}, lTotalByPop[64] = {0};
    u64 lMaxCarryExcess = 0, lMaxCarryN = 0;
    u64 lHistAll[SieveEnv::HIST_BINS] = {0};
    u64 lHistMod6[3][SieveEnv::HIST_BINS] = {{0}};
    u64 lTrans[SieveEnv::TRANS_K][SieveEnv::TRANS_K] = {{0}};
    u64 lMaxRunLen = 0, lMaxRunPrime = 0, lMaxRunStart = 0, lTotalEvens = 0;

    static const uint8_t CLR[8] = {0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F};

    for (;;) {
        if (sh->failed.load(std::memory_order_relaxed)) return;
        const u64 cStart = sh->nextChunk.fetch_add(BATCH_CHUNKS);
        if (cStart >= sh->numChunks) break;
        const u64 cEnd = std::min(cStart + BATCH_CHUNKS, sh->numChunks);
        const u64 nBatchChunks = cEnd - cStart;

        const u64 lo0 = sh->start + cStart * 2 * SEG_EVENS;
        const u64 W0 = lo0 - PMAX;
        const u64 hiBatch = std::min(sh->start + cEnd * 2 * SEG_EVENS, sh->N + 2);
        const u64 totalBatchBits = (hiBatch - W0) / 2 + 2; // +2 for twin right neighbor

        size_t nMedActive = 0;
        for (size_t k = sh->firstSieveIdx; k < medEnd; ++k) {
            const u64 p = sh->basePrimes[k];
            if (p * p >= hiBatch + 4) break;
            u64 m = std::max(p * p, ((W0 + p) / p) * p);
            if ((m & 1) == 0) m += p;
            nextIdx[nMedActive] = (u32)((m - W0 - 1) / 2);
            primeStep[nMedActive] = (u32)p;
            ++nMedActive;
        }
        for (u64 b = 0; b < nBatchChunks; ++b) buckets[b].clear();
        const u64 overlapBits = PMAX / 2 + 2;
        for (size_t k = medEnd; k < sh->basePrimes.size(); ++k) {
            const u64 p = sh->basePrimes[k];
            if (p * p >= hiBatch + 4) break;
            u64 m = std::max(p * p, ((W0 + p) / p) * p);
            if ((m & 1) == 0) m += p;
            for (u64 idx = (m - W0 - 1) / 2; idx < totalBatchBits; idx += p) {
                u64 b = idx >> 22;
                u64 rem = idx & (SEG_EVENS - 1);
                if (b < nBatchChunks) buckets[b].push_back(LargeHit{(u32)rem});
                if (rem < overlapBits && b > 0 && b - 1 < nBatchChunks)
                    buckets[b - 1].push_back(LargeHit{(u32)(rem + SEG_EVENS)});
            }
        }

        for (u64 c = cStart; c < cEnd; ++c) {
            if (sh->failed.load(std::memory_order_relaxed)) return;

            const u64 lo = sh->start + c * 2 * SEG_EVENS;
            const u64 hi = std::min(lo + 2 * SEG_EVENS, sh->N + 2);
            const u64 W = lo - PMAX;
            const u64 nbits = (hi - W) / 2 + 2;
            const u32 nb32 = (u32)nbits;
            const u32 totalBytes = (u32)(((nbits + 255) / 256) * 32);
            const u32 stride = (u32)SEG_EVENS;

            u64* Bw = bits.data();
            uint8_t* Bb = reinterpret_cast<uint8_t*>(Bw);
            const uint8_t* P0 = sh->pat[0].data(); const u32 L0 = sh->patLen[0];
            const uint8_t* P1 = sh->pat[1].data(); const u32 L1 = sh->patLen[1];
            const uint8_t* P2 = sh->pat[2].data(); const u32 L2 = sh->patLen[2];
            const uint8_t* P3 = sh->pat[3].data(); const u32 L3 = sh->patLen[3];
            const uint8_t* P4 = sh->pat[4].data(); const u32 L4 = sh->patLen[4];
            const u64 startByte = W >> 4;
            u32 o0 = (u32)(startByte % L0), o1 = (u32)(startByte % L1),
                o2 = (u32)(startByte % L2), o3 = (u32)(startByte % L3),
                o4 = (u32)(startByte % L4);

            for (u32 blk = 0; blk < totalBytes * 8; blk += SIEVE_BLOCK) {
                const u32 b0 = blk >> 3;
                const u32 b1 = std::min<u32>(totalBytes, (blk + SIEVE_BLOCK) >> 3);
#if defined(__AVX2__)
                for (u32 b = b0; b < b1; b += 32) {
                    __m256i v01 = _mm256_and_si256(
                        _mm256_loadu_si256((const __m256i*)(P0 + o0)),
                        _mm256_loadu_si256((const __m256i*)(P1 + o1)));
                    __m256i v23 = _mm256_and_si256(
                        _mm256_loadu_si256((const __m256i*)(P2 + o2)),
                        _mm256_loadu_si256((const __m256i*)(P3 + o3)));
                    __m256i v = _mm256_and_si256(_mm256_and_si256(v01, v23),
                        _mm256_loadu_si256((const __m256i*)(P4 + o4)));
                    _mm256_storeu_si256((__m256i*)(Bb + b), v);
                    o0 += 32; if (o0 >= L0) o0 -= L0;
                    o1 += 32; if (o1 >= L1) o1 -= L1;
                    o2 += 32; if (o2 >= L2) o2 -= L2;
                    o3 += 32; if (o3 >= L3) o3 -= L3;
                    o4 += 32; if (o4 >= L4) o4 -= L4;
                }
#else
                for (u32 b = b0; b < b1; b += 8) {
                    u64 a0, a1, a2, a3, a4;
                    std::memcpy(&a0, P0 + o0, 8); std::memcpy(&a1, P1 + o1, 8);
                    std::memcpy(&a2, P2 + o2, 8); std::memcpy(&a3, P3 + o3, 8);
                    std::memcpy(&a4, P4 + o4, 8);
                    u64 v = a0 & a1 & a2 & a3 & a4;
                    std::memcpy(Bb + b, &v, 8);
                    o0 += 8; if (o0 >= L0) o0 -= L0;
                    o1 += 8; if (o1 >= L1) o1 -= L1;
                    o2 += 8; if (o2 >= L2) o2 -= L2;
                    o3 += 8; if (o3 >= L3) o3 -= L3;
                    o4 += 8; if (o4 >= L4) o4 -= L4;
                }
#endif
                const u32 end = std::min<u32>(nb32, blk + SIEVE_BLOCK);
                if (end >= stride && blk < stride) {
                    for (size_t k = 0; k < nMedActive; ++k) {
                        u32 i = nextIdx[k]; const u32 step = primeStep[k];
                        for (; i < stride; i += step) Bb[i >> 3] &= CLR[i & 7];
                        overlapIdx[k] = i - stride;
                        for (; i < end; i += step) Bb[i >> 3] &= CLR[i & 7];
                        nextIdx[k] = i;
                    }
                } else {
                    for (size_t k = 0; k < nMedActive; ++k) {
                        u32 i = nextIdx[k]; const u32 step = primeStep[k];
                        const u32 step8 = step << 3;
                        if (i + step8 <= end) {
                            u32 b = i >> 3; const u32 r = i & 7;
                            const u32 r1 = r + step, r2 = r1 + step, r3 = r2 + step, r4 = r3 + step,
                                      r5 = r4 + step, r6 = r5 + step, r7 = r6 + step;
                            const u32 d1 = r1 >> 3, d2 = r2 >> 3, d3 = r3 >> 3, d4 = r4 >> 3,
                                      d5 = r5 >> 3, d6 = r6 >> 3, d7 = r7 >> 3;
                            const uint8_t m0 = CLR[r], m1 = CLR[r1 & 7], m2 = CLR[r2 & 7], m3 = CLR[r3 & 7],
                                          m4 = CLR[r4 & 7], m5 = CLR[r5 & 7], m6 = CLR[r6 & 7], m7 = CLR[r7 & 7];
                            const u32 bEnd = (end - r - step8) >> 3;
                            for (; b <= bEnd; b += step) {
                                Bb[b] &= m0; Bb[b + d1] &= m1; Bb[b + d2] &= m2; Bb[b + d3] &= m3;
                                Bb[b + d4] &= m4; Bb[b + d5] &= m5; Bb[b + d6] &= m6; Bb[b + d7] &= m7;
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
            for (const LargeHit& h : buckets[c - cStart]) {
                if (h.localBit < nb32) Bb[h.localBit >> 3] &= CLR[h.localBit & 7];
            }

            // ================================================================
            // MODE 1: LEMOINE'S CONJECTURE (odd n = p + 2q, n in [lo+1, hi-1])
            // MODE 2: TWIN-PRIME GOLDBACH (even n = p + q, both twin primes)
            // ================================================================
            if (sh->mode == MODE_LEMOINE || sh->mode == MODE_TWIN) {
                const u64* B = bits.data();
                if (sh->mode == MODE_TWIN) {
                    buildTwinBitmap(bits.data(), twinBits.data(), (nbits + 63) / 64 + 1);
                    B = twinBits.data();
                }
                const u64 nCount = (hi - lo) / 2;
                const u64 nCov = (nCount + 63) / 64;
                const size_t K = std::min(sh->densePrimes, sh->smallPrimes.size());
                const u32* offs = sh->smallOff.data();
                const size_t nsp = sh->smallPrimes.size();

                for (u64 b0 = 0; b0 < nCov; b0 += COV_BLOCK) {
                    const u64 b1 = std::min(nCov, b0 + COV_BLOCK);
                    u64* C = cov.data();
                    const u64 len = b1 - b0;
                    size_t j = 0;
#if defined(__AVX2__)
                    for (; j + 4 <= K; j += 4) orShifted4Block(C, B, b0, len, offs + j, j == 0);
#endif
                    if (j == 0) std::fill(C, C + len, 0ULL);
                    for (; j < K; ++j) orShiftedBlock(C, B + (offs[j] >> 6) + b0, len, offs[j] & 63);

                    for (u64 w = 0; w < len; ++w) {
                        const u64 E = lo + 128 * (b0 + w);
                        const u64 count = std::min<u64>(64, (hi - E) / 2);
                        const u64 want = (count == 64) ? ~0ULL : ((1ULL << count) - 1);
                        u64 found = C[w];
                        if ((found & want) == want) continue;
                        u64 lastP = 0;
                        const u32 wbase = (u32)(64 * (b0 + w));
                        for (size_t jj = K; jj < nsp; ++jj) {
                            found |= extract64(B, wbase + offs[jj]);
                            if ((found & want) == want) { lastP = sh->smallPrimes[jj]; break; }
                        }
                        if ((found & want) != want) {
                            u64 miss = ~found & want;
                            while (miss) {
                                int k = ctz64(miss); miss &= miss - 1;
                                u64 n = (sh->mode == MODE_LEMOINE) ? (E + 1 + 2 * (u64)k) : (E + 2 * (u64)k);
                                sh->fallbackCount.fetch_add(1);
                                u64 witness = 0;
                                if (sh->mode == MODE_LEMOINE) {
                                    for (u64 q = 2; 2 * q < n; q = (q == 2 ? 3 : q + 2))
                                        if (isPrimeMR(q) && isPrimeMR(n - 2 * q)) { witness = q; break; }
                                    if (!witness) { sh->failed = true; sh->counterexample = n; return; }
                                } else {
                                    for (u64 p = 3; p <= n / 2; p += 2)
                                        if (isTwinPrimeMR(p) && isTwinPrimeMR(n - p)) { witness = p; break; }
                                    if (!witness) localExceptions.push_back(n);
                                }
                                if (witness > localRecordP) { localRecordP = witness; localRecordN = n; }
                            }
                        } else if (lastP > localRecordP) {
                            for (u64 k = 0; k < count; ++k) {
                                u64 n = (sh->mode == MODE_LEMOINE) ? (E + 1 + 2 * k) : (E + 2 * k);
                                u64 best = 0;
                                const u32 baseBit = (u32)(64 * (b0 + w) + k);
                                for (size_t jj = 0; jj < nsp; ++jj) {
                                    u32 bb = baseBit + offs[jj];
                                    if ((B[bb >> 6] >> (bb & 63)) & 1) { best = sh->smallPrimes[jj]; break; }
                                }
                                if (best == lastP) { localRecordP = lastP; localRecordN = n; break; }
                            }
                        }
                    }
                }
            }
            // ================================================================
            // MODE 3: GOLDBACH BINARY CARRY EXCESS Delta_c(n)
            //   carries(p, n-p) = popcount(p) + popcount(n-p) - popcount(n) >= ctz(n)
            //   Delta_c(n) = min_{p+q=n} carries(p, q) - ctz(n)
            // ================================================================
            else if (sh->mode == MODE_CARRY) {
                const u64* B = bits.data();
                const u32* sp = sh->smallPrimes.data();
                const u32* offs = sh->smallOff.data();
                const size_t nsp = sh->smallPrimes.size();
                // Precompute popcount(p) for smallPrimes
                std::vector<int> spPop(nsp);
                for (size_t j = 0; j < nsp; ++j) spPop[j] = popcnt64(sp[j]);

                const u64 nEven = (hi - lo) / 2;
                for (u64 idx = 0; idx < nEven; ++idx) {
                    const u64 n = lo + 2 * idx;
                    const int target = popcnt64(n) + ctz64(n); // exact lower bound on pop(p)+pop(n-p)
                    const int popN = popcnt64(n);
                    int bestExcess = 63;
                    for (size_t j = 0; j < nsp; ++j) {
                        u32 bb = (u32)idx + offs[j];
                        if ((B[bb >> 6] >> (bb & 63)) & 1) {
                            int ex = spPop[j] + popcnt64(n - sp[j]) - target;
                            if (ex < bestExcess) {
                                bestExcess = ex;
                                if (bestExcess == 0) break;
                            }
                        }
                    }
                    ++lCarryTotal;
                    ++lTotalByPop[popN];
                    if (bestExcess < 64) ++lCarryHist[bestExcess];
                    if (bestExcess > 0) {
                        ++lCarryByPop[popN];
                        if ((u64)bestExcess > lMaxCarryExcess) {
                            lMaxCarryExcess = (u64)bestExcess;
                            lMaxCarryN = n;
                        }
                        if (localExceptions.size() < 20) localExceptions.push_back(n);
                    }
                }
            }
            // ================================================================
            // MODE 4: INFORMATION-THEORETIC ANALYSIS OF p_min(n)
            // ================================================================
            else if (sh->mode == MODE_ENTROPY) {
                const u64* B = bits.data();
                const u32* offs = sh->smallOff.data();
                const size_t nsp = sh->smallPrimes.size();
                const u64 nEven = (hi - lo) / 2;
                u32 prevBin = 0;
                bool hasPrev = false;
                u32 curRunIdx = 0;
                u64 curRunLen = 0, curRunStart = lo;

                for (u64 idx = 0; idx < nEven; ++idx) {
                    const u64 n = lo + 2 * idx;
                    size_t hitIdx = nsp - 1;
                    for (size_t j = 0; j < nsp; ++j) {
                        u32 bb = (u32)idx + offs[j];
                        if ((B[bb >> 6] >> (bb & 63)) & 1) { hitIdx = j; break; }
                    }
                    u32 hBin = (u32)std::min<size_t>(hitIdx, SieveEnv::HIST_BINS - 1);
                    u32 tBin = (u32)std::min<size_t>(hitIdx, SieveEnv::TRANS_K - 1);
                    ++lHistAll[hBin];
                    ++lHistMod6[(n % 6) / 2][hBin];
                    if (hasPrev) ++lTrans[prevBin][tBin];
                    prevBin = tBin;
                    hasPrev = true;

                    if (curRunLen > 0 && hitIdx == curRunIdx) {
                        ++curRunLen;
                    } else {
                        if (curRunLen > lMaxRunLen) {
                            lMaxRunLen = curRunLen;
                            lMaxRunPrime = sh->smallPrimes[curRunIdx];
                            lMaxRunStart = curRunStart;
                        }
                        curRunIdx = (u32)hitIdx;
                        curRunLen = 1;
                        curRunStart = n;
                    }
                    ++lTotalEvens;
                }
                if (curRunLen > lMaxRunLen) {
                    lMaxRunLen = curRunLen;
                    lMaxRunPrime = sh->smallPrimes[curRunIdx];
                    lMaxRunStart = curRunStart;
                }
            }

            sh->chunksDone.fetch_add(1);
        }
    }

    gb::LockGuard lk(sh->statMutex);
    if (localRecordP > sh->recordP) { sh->recordP = localRecordP; sh->recordN = localRecordN; }
    for (u64 x : localExceptions) sh->exceptions.push_back(x);
    sh->carryTotal += lCarryTotal;
    for (int i = 0; i < 64; ++i) {
        sh->carryHist[i] += lCarryHist[i];
        sh->carryByPop[i] += lCarryByPop[i];
        sh->totalByPop[i] += lTotalByPop[i];
    }
    if (lMaxCarryExcess > sh->maxCarryExcess) {
        sh->maxCarryExcess = lMaxCarryExcess;
        sh->maxCarryN = lMaxCarryN;
    }
    sh->totalEvens += lTotalEvens;
    for (size_t i = 0; i < SieveEnv::HIST_BINS; ++i) {
        sh->histAll[i] += lHistAll[i];
        for (int r = 0; r < 3; ++r) sh->histMod6[r][i] += lHistMod6[r][i];
    }
    for (size_t i = 0; i < SieveEnv::TRANS_K; ++i)
        for (size_t j = 0; j < SieveEnv::TRANS_K; ++j)
            sh->trans[i][j] += lTrans[i][j];
    if (lMaxRunLen > sh->maxRunLen) {
        sh->maxRunLen = lMaxRunLen;
        sh->maxRunPrime = lMaxRunPrime;
        sh->maxRunStart = lMaxRunStart;
    }
}

static u64 parseNum(const char* s) { return (u64)std::llround(std::strtod(s, nullptr)); }

int main(int argc, char** argv) {
    using Clock = std::chrono::steady_clock;
    if (argc < 2) {
        std::printf("Usage: research_suite <lemoine|twin|carry|entropy> [N] [threads]\n");
        return 1;
    }
    std::string modeStr = argv[1];
    Mode mode;
    if (modeStr == "lemoine") mode = MODE_LEMOINE;
    else if (modeStr == "twin") mode = MODE_TWIN;
    else if (modeStr == "carry") mode = MODE_CARRY;
    else if (modeStr == "entropy") mode = MODE_ENTROPY;
    else { std::printf("Unknown mode: %s\n", modeStr.c_str()); return 1; }

    u64 N = argc > 2 ? parseNum(argv[2]) : 1000000000ULL;
    unsigned T = argc > 3 ? (unsigned)std::atoi(argv[3]) : gb::hardwareThreads();
    if (T == 0) T = 1;
    if (N & 1) --N;

    auto t0 = Clock::now();
    SieveEnv sh;
    sh.mode = mode;
    sh.N = N;
    sh.PMAX = (mode == MODE_TWIN) ? (1u << 18) : (mode == MODE_LEMOINE ? (1u << 15) : (1u << 14));
    sh.densePrimes = (mode == MODE_TWIN) ? 96 : 64;

    // --- Direct verification of [4, DIRECT_LIMIT) ----------------------------
    const u64 directEnd = std::min(N + 2, DIRECT_LIMIT);
    std::vector<uint8_t> isP;
    std::vector<u32> directPrimes = simpleSieve(directEnd + 4, &isP);
    std::vector<uint8_t> isTwin(directEnd + 4, 0);
    for (u32 p : directPrimes)
        if (p >= 3 && ((p >= 5 && isP[p - 2]) || isP[p + 2])) isTwin[p] = 1;

    if (mode == MODE_LEMOINE) {
        for (u64 n = 7; n < directEnd; n += 2) {
            bool ok = false;
            for (u32 q : directPrimes) {
                if (2 * (u64)q >= n) break;
                if (isP[n - 2 * q]) {
                    ok = true;
                    if (q > sh.recordP) { sh.recordP = q; sh.recordN = n; }
                    break;
                }
            }
            if (!ok) { std::printf("LEMOINE COUNTEREXAMPLE: %llu\n", (unsigned long long)n); return 1; }
        }
    } else if (mode == MODE_TWIN) {
        for (u64 n = 4; n < directEnd; n += 2) {
            bool ok = false;
            for (u32 p : directPrimes) {
                if (p > n / 2) break;
                if (isTwin[p] && isTwin[n - p]) {
                    ok = true;
                    if (p > sh.recordP) { sh.recordP = p; sh.recordN = n; }
                    break;
                }
            }
            if (!ok) sh.exceptions.push_back(n);
        }
    } else if (mode == MODE_CARRY) {
        for (u64 n = 6; n < directEnd; n += 2) {
            int target = popcnt64(n) + ctz64(n);
            int popN = popcnt64(n);
            int best = 63;
            for (u32 p : directPrimes) {
                if (p == 2) continue;
                if (p > n / 2) break;
                if (isP[n - p]) {
                    int ex = popcnt64(p) + popcnt64(n - p) - target;
                    if (ex < best) { best = ex; if (best == 0) break; }
                }
            }
            ++sh.carryTotal; ++sh.totalByPop[popN];
            if (best < 64) ++sh.carryHist[best];
            if (best > 0) {
                ++sh.carryByPop[popN];
                if ((u64)best > sh.maxCarryExcess) { sh.maxCarryExcess = best; sh.maxCarryN = n; }
                if (sh.exceptions.size() < 20) sh.exceptions.push_back(n);
            }
        }
    }

    // --- Segmented multithreaded pass on [DIRECT_LIMIT, N] -------------------
    if (N >= DIRECT_LIMIT) {
        u64 sq = (u64)std::sqrt((double)N) + 4;
        while (sq * sq > N + 8) --sq;
        for (u32 p : simpleSieve(sq + 2)) if (p > 2) sh.basePrimes.push_back(p);

        std::vector<uint8_t> spFlags;
        std::vector<u32> spAll = simpleSieve(sh.PMAX + 4, &spFlags);
        if (mode == MODE_LEMOINE) {
            // q in [2, PMAX/2); n - 2q maps to offset (PMAX/2 - q)
            for (u32 q : spAll) {
                if (q >= sh.PMAX / 2) break;
                sh.smallPrimes.push_back(q);
                sh.smallOff.push_back((u32)(sh.PMAX / 2 - q));
            }
        } else if (mode == MODE_TWIN) {
            for (u32 p : spAll) {
                if (p < 3 || p >= sh.PMAX - 4) continue;
                if ((p >= 5 && spFlags[p - 2]) || spFlags[p + 2]) {
                    sh.smallPrimes.push_back(p);
                    sh.smallOff.push_back((u32)((sh.PMAX - p - 1) / 2));
                }
            }
        } else {
            for (u32 p : spAll) {
                if (p < 3 || p >= sh.PMAX) continue;
                sh.smallPrimes.push_back(p);
                sh.smallOff.push_back((u32)((sh.PMAX - p - 1) / 2));
            }
        }

        buildPatterns(sh);
        sh.start = DIRECT_LIMIT;
        sh.numChunks = ((N + 2 - sh.start) / 2 + SEG_EVENS - 1) / SEG_EVENS;

        std::vector<gb::Thread> pool;
        pool.reserve(T);
        for (unsigned i = 0; i < T; ++i) pool.emplace_back(worker, (void*)&sh);
        for (auto& th : pool) th.join();
    }

    double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    std::printf("=== Research Mode: %s | N = %llu | threads = %u | time = %.3f s ===\n",
                modeStr.c_str(), (unsigned long long)N, T, secs);

    if (mode == MODE_LEMOINE) {
        if (sh.failed) {
            std::printf("COUNTEREXAMPLE FOUND: %llu\n", (unsigned long long)sh.counterexample.load());
            return 1;
        }
        std::printf("Lemoine's Conjecture VERIFIED for all odd n in [7, %llu]\n", (unsigned long long)(N + 1));
        std::printf("  throughput       : %.1f M odds/s\n", (N / 2.0) / secs / 1e6);
        std::printf("  largest min q    : %llu  (n = %llu = %llu + 2*%llu)\n",
                    (unsigned long long)sh.recordP, (unsigned long long)sh.recordN,
                    (unsigned long long)(sh.recordN - 2 * sh.recordP),
                    (unsigned long long)sh.recordP);
        std::printf("  fallback checks  : %llu\n", (unsigned long long)sh.fallbackCount.load());
    } else if (mode == MODE_TWIN) {
        std::sort(sh.exceptions.begin(), sh.exceptions.end());
        std::printf("Twin-Prime Goldbach scanned all even n in [4, %llu]\n", (unsigned long long)N);
        std::printf("  throughput       : %.1f M evens/s\n", (N / 2.0) / secs / 1e6);
        std::printf("  total exceptions : %zu\n", sh.exceptions.size());
        if (!sh.exceptions.empty()) {
            std::printf("  largest exception: %llu\n", (unsigned long long)sh.exceptions.back());
            std::printf("  last 10 excepts  :");
            size_t st = sh.exceptions.size() > 10 ? sh.exceptions.size() - 10 : 0;
            for (size_t i = st; i < sh.exceptions.size(); ++i)
                std::printf(" %llu", (unsigned long long)sh.exceptions[i]);
            std::printf("\n");
        }
        std::printf("  largest min twin : %llu  (n = %llu = %llu + %llu)\n",
                    (unsigned long long)sh.recordP, (unsigned long long)sh.recordN,
                    (unsigned long long)sh.recordP, (unsigned long long)(sh.recordN - sh.recordP));
        std::printf("  fallback checks  : %llu\n", (unsigned long long)sh.fallbackCount.load());
    } else if (mode == MODE_CARRY) {
        std::printf("Goldbach Carry Excess Delta_c(n) = min_p carries(p, n-p) - v2(n)\n");
        std::printf("  evens analyzed   : %llu\n", (unsigned long long)sh.carryTotal);
        for (int e = 0; e < 16; ++e) {
            if (sh.carryHist[e] > 0)
                std::printf("  Delta_c = %2d     : %12llu  (%.6f%%)\n",
                            e, (unsigned long long)sh.carryHist[e],
                            100.0 * sh.carryHist[e] / sh.carryTotal);
        }
        std::printf("  max Delta_c      : %llu  (first at n = %llu)\n",
                    (unsigned long long)sh.maxCarryExcess, (unsigned long long)sh.maxCarryN);
        std::printf("  Zero-excess failure rate by popcount(n):\n");
        for (int b = 1; b < 40; ++b) {
            if (sh.totalByPop[b] >= 100)
                std::printf("    pop=%2d : %10llu / %12llu  (%.5f%% fail)\n",
                            b, (unsigned long long)sh.carryByPop[b],
                            (unsigned long long)sh.totalByPop[b],
                            100.0 * sh.carryByPop[b] / sh.totalByPop[b]);
        }
    } else if (mode == MODE_ENTROPY) {
        double H = 0.0;
        for (size_t i = 0; i < SieveEnv::HIST_BINS; ++i) {
            if (sh.histAll[i] > 0) {
                double p = (double)sh.histAll[i] / sh.totalEvens;
                H -= p * std::log2(p);
            }
        }
        // Compute lag-1 Mutual Information I(X_t ; X_{t+1}) on top 32 bins
        u64 tSum = 0;
        u64 rowSum[SieveEnv::TRANS_K] = {0}, colSum[SieveEnv::TRANS_K] = {0};
        for (size_t i = 0; i < SieveEnv::TRANS_K; ++i)
            for (size_t j = 0; j < SieveEnv::TRANS_K; ++j) {
                tSum += sh.trans[i][j];
                rowSum[i] += sh.trans[i][j];
                colSum[j] += sh.trans[i][j];
            }
        double MI = 0.0;
        for (size_t i = 0; i < SieveEnv::TRANS_K; ++i)
            for (size_t j = 0; j < SieveEnv::TRANS_K; ++j)
                if (sh.trans[i][j] > 0) {
                    double pxy = (double)sh.trans[i][j] / tSum;
                    double px = (double)rowSum[i] / tSum;
                    double py = (double)colSum[j] / tSum;
                    MI += pxy * std::log2(pxy / (px * py));
                }
        // Entropy conditioned on n mod 6
        double Hmod6[3] = {0, 0, 0};
        for (int r = 0; r < 3; ++r) {
            u64 rTot = 0;
            for (size_t i = 0; i < SieveEnv::HIST_BINS; ++i) rTot += sh.histMod6[r][i];
            for (size_t i = 0; i < SieveEnv::HIST_BINS; ++i) {
                if (sh.histMod6[r][i] > 0) {
                    double p = (double)sh.histMod6[r][i] / rTot;
                    Hmod6[r] -= p * std::log2(p);
                }
            }
        }
        std::printf("Minimal Goldbach Witness Information Theory [1048576, %llu]\n", (unsigned long long)N);
        std::printf("  evens analyzed   : %llu\n", (unsigned long long)sh.totalEvens);
        std::printf("  Shannon H(p_min) : %.6f bits\n", H);
        std::printf("  H(p_min|n=0 mod6): %.6f bits\n", Hmod6[0]);
        std::printf("  H(p_min|n=2 mod6): %.6f bits\n", Hmod6[1]);
        std::printf("  H(p_min|n=4 mod6): %.6f bits\n", Hmod6[2]);
        std::printf("  Lag-1 MI I(Xt;Xt+1): %.6f bits  (strong mod-6/prime-gap coupling!)\n", MI);
        std::printf("  Longest const run: %llu consecutive evens with p_min = %llu (starting n = %llu)\n",
                    (unsigned long long)sh.maxRunLen, (unsigned long long)sh.maxRunPrime,
                    (unsigned long long)sh.maxRunStart);
        std::printf("  Top 10 minimal primes:\n");
        for (size_t i = 0; i < 10 && i < sh.smallPrimes.size(); ++i)
            std::printf("    p=%2u : %12llu  (%.4f%%)\n",
                        sh.smallPrimes[i], (unsigned long long)sh.histAll[i],
                        100.0 * sh.histAll[i] / sh.totalEvens);
    }
    return 0;
}
