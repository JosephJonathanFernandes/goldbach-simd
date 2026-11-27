// gosper_hunt.cpp - Combinatorial Gosper-Hack Hunter for Minimal Carry-Optimal
//                   Goldbach Exceptions n_min(w) (w = 1, 2, ..., 13+)
//
// Defines the new integer sequence:
//   a(w) = n_min(w) = smallest even integer n >= 6 with popcount(n - 2) = w
//                     such that Delta_c*(n) > 0 (no carry-optimal Goldbach pair).
//
// Uses bit-shell decomposition (highest bit B of M = n - 2) + Gosper's hack to
// visit ONLY integers of exact binary weight w, O(1) odd-prime bitset lookup
// below 2^27, wheel trial division, and fast x87 80-bit mulmod + Sinclair 7-base
// deterministic 64-bit Miller-Rabin.
//
// Build: g++ -O3 -march=native -std=c++14 -Wall -Wno-format gosper_hunt.cpp -o gosper_hunt
// Usage: gosper_hunt [max_w] [max_bits] [threads]

#define __USE_MINGW_ANSI_STDIO 1
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

typedef uint64_t u64;
typedef uint32_t u32;

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
#else
typedef std::mutex Mutex;
class Thread {
    std::thread t_;
public:
    Thread(void (*fn)(void*), void* arg) : t_(fn, arg) {}
    void join() { t_.join(); }
};
inline unsigned hardwareThreads() { return std::thread::hardware_concurrency(); }
#endif
struct LockGuard {
    Mutex& m;
    explicit LockGuard(Mutex& mm) : m(mm) { m.lock(); }
    ~LockGuard() { m.unlock(); }
};
}  // namespace gb

static inline int ctz64(u64 x) { return __builtin_ctzll(x); }

// Reference slow mulmod (used for self-test & fallback >= 2^62)
static u64 mulmodSlow(u64 a, u64 b, u64 m) {
    u64 r = 0; a %= m;
    while (b) {
        if (b & 1) { r += a; if (r >= m || r < a) r -= m; }
        u64 a2 = a + a; if (a2 >= m || a2 < a) a2 -= m; a = a2;
        b >>= 1;
    }
    return r;
}

// Fast mulmod for m < 2^62 using x87 80-bit long double (64-bit significand) or __int128
static inline u64 mulmodFast(u64 a, u64 b, u64 m) {
#if defined(__SIZEOF_INT128__)
    return (u64)((unsigned __int128)a * b % m);
#else
    if (sizeof(long double) >= 10 && m < (1ULL << 62)) {
        u64 q = (u64)(((long double)a * (long double)b) / (long double)m);
        int64_t r = (int64_t)(a * b - q * m);
        while (r < 0) r += (int64_t)m;
        while ((u64)r >= m) r -= (int64_t)m;
        return (u64)r;
    }
    return mulmodSlow(a, b, m);
#endif
}

static inline u64 powmodFast(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmodFast(r, a, m); a = mulmodFast(a, a, m); e >>= 1; }
    return r;
}

static inline bool mrWitness(u64 a, u64 d, int s, u64 n) {
    a %= n;
    if (a == 0) return true;
    u64 x = powmodFast(a, d, n);
    if (x == 1 || x == n - 1) return true;
    for (int r = 1; r < s; ++r) {
        x = mulmodFast(x, x, n);
        if (x == n - 1) return true;
    }
    return false;
}

static const u32 BITSET_MAX = 1u << 27; // 134,217,728 (8 MB odd-only bitset)

struct PrimeOracle {
    std::vector<u64> oddBits;

    void init() {
        const u32 nOdds = BITSET_MAX / 2 + 1;
        oddBits.assign(nOdds / 64 + 2, ~0ULL);
        oddBits[0] &= ~1ULL; // 1 is not prime
        u64* B = oddBits.data();
        const u32 sq = (u32)std::sqrt((double)BITSET_MAX) + 1;
        for (u32 p = 3; p <= sq; p += 2) {
            if ((B[(p >> 1) >> 6] >> ((p >> 1) & 63)) & 1ULL) {
                for (u32 m = p * p; m <= BITSET_MAX; m += 2 * p) {
                    u32 i = m >> 1;
                    B[i >> 6] &= ~(1ULL << (i & 63));
                }
            }
        }
    }

    inline bool isOddPrime(u64 x) const {
        if (x < BITSET_MAX) {
            u32 i = (u32)(x >> 1);
            return (oddBits[i >> 6] >> (i & 63)) & 1ULL;
        }
        // Quick trial division by small odd primes (~83% filtered before Miller-Rabin)
        if (x % 3 == 0 || x % 5 == 0 || x % 7 == 0 || x % 11 == 0 || x % 13 == 0) return false;
        if (x % 17 == 0 || x % 19 == 0 || x % 23 == 0 || x % 29 == 0 || x % 31 == 0) return false;
        if (x % 37 == 0 || x % 41 == 0 || x % 43 == 0 || x % 47 == 0 || x % 53 == 0) return false;

        u64 d = x - 1;
        int s = ctz64(d);
        d >>= s;
        // Deterministic for all n < 4,759,123,141 (Jaeschke): just 3 bases {2, 7, 61}!
        if (x < 4759123141ULL) {
            return mrWitness(2, d, s, x) && mrWitness(7, d, s, x) && mrWitness(61, d, s, x);
        }
        // Deterministic for all 64-bit integers (Sinclair 2011): 7 bases
        static const u64 sinclair[7] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
        for (int i = 0; i < 7; ++i)
            if (!mrWitness(sinclair[i], d, s, x)) return false;
        return true;
    }

    // Exact check: does n = M + 2 (where M has highest bit `hiBit` and lower bits `M_lo`)
    // admit ANY carry-optimal Goldbach pair (1|s, 1|(M^s))?
    inline bool hasZeroExcess(u64 M, u64 M_lo) const {
        if (M_lo == 0) return false;
        for (u64 s = M_lo & (0ULL - M_lo); s != 0; s = (s - M_lo) & M_lo) {
            if (isOddPrime(1ULL | s) && isOddPrime(1ULL | (M ^ s)))
                return true;
        }
        return false;
    }
};

// Advance `x` to the next integer with the exact same popcount (Gosper's hack).
static inline u64 gosperNext(u64 x) {
    u64 c = x & (0ULL - x);
    u64 r = x + c;
    return (((r ^ x) >> 2) / c) | r;
}

// Binomial coefficient C(n, k)
static u64 binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    u64 r = 1;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

// Unrank the `idx`-th (0-indexed) k-combination of {0, 1, ..., n-1} in colexicographic
// (= numeric increasing!) order! This lets us partition Gosper's hack cleanly across T threads!
static u64 unrankColex(int n, int k, u64 idx) {
    u64 mask = 0;
    int remK = k;
    for (int bit = n - 1; bit >= 0 && remK > 0; --bit) {
        u64 c = binom(bit, remK);
        if (idx >= c) {
            idx -= c;
            mask |= (1ULL << bit);
            --remK;
        }
    }
    return mask;
}

struct ShellJob {
    const PrimeOracle* oracle;
    int B;              // highest bit position of M (in 1..60); M = (1ULL << B) | (comb << 1)
    int w;              // total weight of M; comb chooses w-1 bits from {0..B-2} (shifted << 1)
    u64 totalCombs;
    std::atomic<u64> nextShard{0};
    u64 shardSize;
    std::atomic<u64> bestM{~0ULL};
};

static void shellWorker(void* arg) {
    ShellJob* job = static_cast<ShellJob*>(arg);
    const u64 hiBit = 1ULL << job->B;
    const int nSlots = job->B - 1;
    const int kChoose = job->w - 1;

    for (;;) {
        u64 s0 = job->nextShard.fetch_add(job->shardSize);
        if (s0 >= job->totalCombs) break;
        // If an earlier shard already found an exception smaller than our shard's minimum, skip!
        u64 comb = unrankColex(nSlots, kChoose, s0);
        if ((hiBit | (comb << 1)) >= job->bestM.load(std::memory_order_relaxed)) break;

        u64 s1 = std::min(job->totalCombs, s0 + job->shardSize);
        for (u64 i = s0; i < s1; ++i) {
            u64 M_lo = comb << 1;
            u64 M = hiBit | M_lo;
            if (M >= job->bestM.load(std::memory_order_relaxed)) break;
            if (!job->oracle->hasZeroExcess(M, M_lo)) {
                // Atomically update minimum M
                u64 prev = job->bestM.load();
                while (M < prev && !job->bestM.compare_exchange_weak(prev, M)) {}
                break; // comb is increasing within this shard, so no smaller M here!
            }
            if (i + 1 < s1) comb = gosperNext(comb);
        }
    }
}

int main(int argc, char** argv) {
    using Clock = std::chrono::steady_clock;
    int maxW = argc > 1 ? std::atoi(argv[1]) : 13;
    int maxBits = argc > 2 ? std::atoi(argv[2]) : 36;
    unsigned T = argc > 3 ? (unsigned)std::atoi(argv[3]) : gb::hardwareThreads();
    if (T == 0) T = 1;

    // 1. Self-test mulmodFast against mulmodSlow
    u64 rng = 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 200000; ++i) {
        rng ^= rng >> 12; rng ^= rng << 25; rng ^= rng >> 27;
        u64 m = ((rng * 0x2545F4914F6CDD1DULL) & ((1ULL << 61) - 1)) | 3ULL;
        u64 a = (rng * 6364136223846793005ULL) % m;
        u64 b = (rng * 1442695040888963407ULL) % m;
        if (mulmodFast(a, b, m) != mulmodSlow(a, b, m)) {
            std::printf("SELF-TEST FAILED on mulmodFast!\n");
            return 2;
        }
    }

    // 1b. Self-test unrankColex against gosperNext
    u64 cTest = (1ULL << 5) - 1;
    for (u64 idx = 0; idx < binom(12, 5); ++idx) {
        if (unrankColex(12, 5, idx) != cTest) {
            std::printf("SELF-TEST FAILED on unrankColex at idx=%llu!\n", (unsigned long long)idx);
            return 2;
        }
        if (idx + 1 < binom(12, 5)) cTest = gosperNext(cTest);
    }

    auto t0 = Clock::now();
    PrimeOracle oracle;
    oracle.init();
    std::printf("PrimeOracle initialized (bitset < 2^27 = %u) in %.3f s | threads = %u\n\n",
                BITSET_MAX, std::chrono::duration<double>(Clock::now() - t0).count(), T);

    std::printf("=== Minimal Carry-Optimal Goldbach Exception Sequence a(w) = n_min(w) ===\n");
    std::printf("  w |            n_min(w) |              hex | bits |   time (s)\n");
    std::printf("----+---------------------+------------------+------+-----------\n");

    for (int w = 1; w <= maxW; ++w) {
        auto tw = Clock::now();
        u64 foundN = 0;
        for (int B = (w == 1 ? 2 : w); B <= maxBits; ++B) {  // n = M + 2 >= 6
            ShellJob job;
            job.oracle = &oracle;
            job.B = B;
            job.w = w;
            job.totalCombs = binom(B - 1, w - 1);
            job.shardSize = std::max<u64>(4096, job.totalCombs / (T * 64));

            if (job.totalCombs < 32768 || T == 1) {
                shellWorker(&job);
            } else {
                std::vector<gb::Thread> pool;
                pool.reserve(T);
                for (unsigned t = 0; t < T; ++t) pool.emplace_back(shellWorker, (void*)&job);
                for (auto& th : pool) th.join();
            }

            u64 bestM = job.bestM.load();
            if (bestM != ~0ULL) {
                foundN = bestM + 2;
                break;
            }
            if (w >= 13 && B >= 30) {
                double el = std::chrono::duration<double>(Clock::now() - tw).count();
                std::fprintf(stderr, "  [w=%d] shell B=%d (%llu candidates) clear in %.2f s\n",
                             w, B, (unsigned long long)job.totalCombs, el);
            }
        }
        double sw = std::chrono::duration<double>(Clock::now() - tw).count();
        if (foundN) {
            int bits = 64 - __builtin_clzll(foundN);
            std::printf(" %2d | %19llu | %16llX |  %2d  | %9.3f\n",
                        w, (unsigned long long)foundN, (unsigned long long)foundN, bits, sw);
        } else {
            std::printf(" %2d | > 2^%-15d |                - | >%2d  | %9.3f\n",
                        w, maxBits + 1, maxBits + 1, sw);
        }
        std::fflush(stdout);
    }
    return 0;
}
