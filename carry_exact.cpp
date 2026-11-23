// carry_exact.cpp - Exact Global Carry-Optimal Goldbach Solver
//
// Theorem: For any even n >= 6 with k = v2(n) = ctz(n) and odd primes p + q = n:
//   carries(p, q) = popcount(p) + popcount(q) - popcount(n) >= k.
//   Equality (Delta_c*(n) = 0, "carry-optimal Goldbach") holds iff
//     p = 1 | s,  q = 1 | ((n - 2) ^ s)
//   for some binary submask s of M = n - 2 with 0 < s < M.
//
// This program computes the EXACT global Carry Excess Delta_c*(n) over ALL primes
// p <= n/2 for every even n in [6, N], using exact submask enumeration for Delta_c* = 0
// and exact prime search when Delta_c*(n) > 0.
//
// Build: g++ -O3 -march=native -std=c++14 -Wall -Wno-format carry_exact.cpp -o carry_exact
// Usage: carry_exact [N] [threads]    (default N = 1e8)

#define __USE_MINGW_ANSI_STDIO 1
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
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

static inline int pop32(u32 x) { return __builtin_popcount(x); }
static inline int ctz32(u32 x) { return __builtin_ctz(x); }
static inline int clz32(u32 x) { return __builtin_clz(x); }

struct Shared {
    u32 N;
    bool zeroOnly = false;
    std::vector<u64> oddPrime;   // bit (x >> 1) is 1 iff odd x is prime (x >= 3; bit 0 = 1 is 0)
    std::vector<u32> primes;     // all odd primes <= N/2 (only populated when !zeroOnly)
    std::vector<uint8_t> pPop;   // popcount(p) for primes[]
    std::atomic<u32> nextLo{6};

    gb::Mutex mtx;
    u64 totalEvens = 0;
    u64 histExcess[32] = {0};
    u64 failByW[32] = {0};       // Delta_c*(n) > 0 grouped by w = popcount(n - 2)
    u64 totalByW[32] = {0};
    u32 maxFailW = 0, maxFailWN = 0;
    u32 maxExcess = 0, maxExcessN = 0;
    u32 lastFailW[32] = {0};     // largest n <= N with Delta_c*(n) > 0 for each w
    std::vector<u32> highWEx;    // exceptions with w >= 11
    std::vector<u32> exGe2;      // numbers with Delta_c*(n) >= 2
};

static inline bool isOddPrime(const u64* B, u32 x) {
    u32 i = x >> 1;
    return (B[i >> 6] >> (i & 63)) & 1ULL;
}

// Exact check: does n admit ANY odd prime pair (p, q) with carries(p, q) == v2(n)?
// Uses the exact theorem: p = 1 | s, q = 1 | (M ^ s) for non-zero submask s of M_lo.
static inline bool hasZeroExcess(const u64* B, u32 n) {
    const u32 M = n - 2;
    if (M == 0) return false;
    const u32 hiBit = 1u << (31 - clz32(M));
    const u32 M_lo = M ^ hiBit;
    if (M_lo == 0) return false; // w = popcount(n - 2) == 1: 0 valid pairs!
    // Enumerate non-zero submasks of M_lo in increasing order (fastest hits!)
    for (u32 s = M_lo & (0u - M_lo); s != 0; s = (s - M_lo) & M_lo) {
        if (isOddPrime(B, 1u | s) && isOddPrime(B, 1u | (M ^ s)))
            return true;
    }
    return false;
}

// Exact global Delta_c*(n) when we already know Delta_c*(n) >= 1.
static int exactGlobalExcess(const Shared* sh, u32 n) {
    const u64* B = sh->oddPrime.data();
    const u32* P = sh->primes.data();
    const uint8_t* PP = sh->pPop.data();
    const size_t np = sh->primes.size();
    const u32 half = n >> 1;
    const int target = pop32(n) + ctz32(n);
    int best = 31;
    for (size_t i = 0; i < np; ++i) {
        const u32 p = P[i];
        if (p > half) break;
        if (isOddPrime(B, n - p)) {
            int ex = (int)PP[i] + pop32(n - p) - target;
            if (ex < best) {
                best = ex;
                if (best == 1) return 1; // cannot be 0 (already checked)
            }
        }
    }
    return best;
}

static void worker(void* arg) {
    Shared* sh = static_cast<Shared*>(arg);
    const u64* B = sh->oddPrime.data();
    const u32 CHUNK = 1u << 17; // 131,072 numbers per chunk

    u64 lTotal = 0, lHist[32] = {0}, lFailW[32] = {0}, lTotW[32] = {0};
    u32 lMaxFailW = 0, lMaxFailWN = 0, lMaxEx = 0, lMaxExN = 0;
    u32 lLastFailW[32] = {0};
    std::vector<u32> lHighW, lExGe2;

    for (;;) {
        u32 lo = sh->nextLo.fetch_add(CHUNK);
        if (lo > sh->N) break;
        u32 hi = std::min<u32>(lo + CHUNK - 2, sh->N);
        for (u32 n = lo; n <= hi; n += 2) {
            int w = pop32(n - 2);
            ++lTotal;
            ++lTotW[w];
            if (hasZeroExcess(B, n)) {
                ++lHist[0];
            } else {
                int ex = sh->zeroOnly ? 1 : exactGlobalExcess(sh, n);
                if (ex < 32) ++lHist[ex];
                ++lFailW[w];
                if (n > lLastFailW[w]) lLastFailW[w] = n;
                if ((u32)w > lMaxFailW) { lMaxFailW = (u32)w; lMaxFailWN = n; }
                if ((u32)ex > lMaxEx) { lMaxEx = (u32)ex; lMaxExN = n; }
                if (w >= 11 && lHighW.size() < 200) lHighW.push_back(n);
                if (!sh->zeroOnly && ex >= 2 && lExGe2.size() < 50) lExGe2.push_back(n);
            }
        }
    }

    gb::LockGuard lk(sh->mtx);
    sh->totalEvens += lTotal;
    for (int i = 0; i < 32; ++i) {
        sh->histExcess[i] += lHist[i];
        sh->failByW[i] += lFailW[i];
        sh->totalByW[i] += lTotW[i];
        if (lLastFailW[i] > sh->lastFailW[i]) sh->lastFailW[i] = lLastFailW[i];
    }
    if (lMaxFailW > sh->maxFailW || (lMaxFailW == sh->maxFailW && lMaxFailWN < sh->maxFailWN)) {
        sh->maxFailW = lMaxFailW; sh->maxFailWN = lMaxFailWN;
    }
    if (lMaxEx > sh->maxExcess || (lMaxEx == sh->maxExcess && lMaxExN < sh->maxExcessN)) {
        sh->maxExcess = lMaxEx; sh->maxExcessN = lMaxExN;
    }
    for (u32 x : lHighW) sh->highWEx.push_back(x);
    for (u32 x : lExGe2) sh->exGe2.push_back(x);
}

int main(int argc, char** argv) {
    using Clock = std::chrono::steady_clock;
    u32 N = argc > 1 ? (u32)std::llround(std::strtod(argv[1], nullptr)) : 100000000u;
    unsigned T = argc > 2 ? (unsigned)std::atoi(argv[2]) : gb::hardwareThreads();
    if (T == 0) T = 1;
    if (N & 1) --N;

    auto t0 = Clock::now();
    Shared sh;
    sh.N = N;
    sh.zeroOnly = (argc > 3 && std::string(argv[3]) == "zero");

    // Sieve odd primes up to N
    const u32 nOdds = N / 2 + 1;
    sh.oddPrime.assign(nOdds / 64 + 2, ~0ULL);
    sh.oddPrime[0] &= ~1ULL; // 1 is not prime
    u64* B = sh.oddPrime.data();
    const u32 sq = (u32)std::sqrt((double)N) + 1;
    for (u32 p = 3; p <= sq; p += 2) {
        if ((B[(p >> 1) >> 6] >> ((p >> 1) & 63)) & 1ULL) {
            for (u32 m = p * p; m <= N; m += 2 * p) {
                u32 i = m >> 1;
                B[i >> 6] &= ~(1ULL << (i & 63));
            }
        }
    }
    if (!sh.zeroOnly) {
        for (u32 p = 3; p <= N / 2; p += 2) {
            if (isOddPrime(B, p)) {
                sh.primes.push_back(p);
                sh.pPop.push_back((uint8_t)pop32(p));
            }
        }
    }

    auto t1 = Clock::now();
    std::vector<gb::Thread> pool;
    pool.reserve(T);
    for (unsigned i = 0; i < T; ++i) pool.emplace_back(worker, (void*)&sh);
    for (auto& th : pool) th.join();

    double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    double sieveSecs = std::chrono::duration<double>(t1 - t0).count();

    std::printf("=== EXACT Global Carry-Optimal Goldbach | N = %u | threads = %u ===\n", N, T);
    std::printf("  time             : %.3f s (sieve %.3f s, exact search %.3f s)\n",
                secs, sieveSecs, secs - sieveSecs);
    std::printf("  evens analyzed   : %llu  ([6, %u])\n", (unsigned long long)sh.totalEvens, N);
    std::printf("  Exact global Delta_c*(n) distribution:\n");
    for (int e = 0; e < 16; ++e) {
        if (sh.histExcess[e] > 0)
            std::printf("    Delta_c* = %2d  : %12llu  (%.6f%%)\n",
                        e, (unsigned long long)sh.histExcess[e],
                        100.0 * sh.histExcess[e] / sh.totalEvens);
    }
    std::printf("  max global Delta_c*: %u  (at n = %u = 0x%X)\n",
                sh.maxExcess, sh.maxExcessN, sh.maxExcessN);
    std::printf("  max w = pop(n-2) with Delta_c* > 0 : w = %u  (at n = %u = 0x%X)\n",
                sh.maxFailW, sh.maxFailWN, sh.maxFailWN);
    std::printf("  Exact exceptions by binary weight w = popcount(n - 2):\n");
    for (int w = 1; w < 31; ++w) {
        if (sh.totalByW[w] > 0)
            std::printf("    w=%2d : %9llu / %11llu  (%9.5f%% fail)  last_n=%u\n",
                        w, (unsigned long long)sh.failByW[w],
                        (unsigned long long)sh.totalByW[w],
                        100.0 * sh.failByW[w] / sh.totalByW[w],
                        sh.lastFailW[w]);
    }
    std::sort(sh.highWEx.begin(), sh.highWEx.end());
    std::printf("  Exceptions with w >= 11 (count=%zu, showing up to 30): ", sh.highWEx.size());
    for (size_t i = 0; i < sh.highWEx.size() && i < 30; ++i)
        std::printf("%u(w=%d) ", sh.highWEx[i], pop32(sh.highWEx[i] - 2));
    std::printf("\n");
    if (!sh.zeroOnly) {
        std::sort(sh.exGe2.begin(), sh.exGe2.end());
        std::printf("  Numbers with Delta_c*(n) >= 2 (showing up to 20): ");
        for (size_t i = 0; i < sh.exGe2.size() && i < 20; ++i) std::printf("%u ", sh.exGe2[i]);
        std::printf("\n");
    }
    return 0;
}
