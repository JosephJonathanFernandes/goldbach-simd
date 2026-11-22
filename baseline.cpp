// baseline.cpp - the straightforward version (full sieve + vector<bool>), timed.
// Build: g++ -O3 -std=c++14 baseline.cpp -o baseline
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

int main(int argc, char** argv) {
    const long long N = argc > 1 ? (long long)atof(argv[1]) : 1000000;
    auto t0 = chrono::steady_clock::now();
    vector<bool> isPrime(N + 1, true);
    isPrime[0] = isPrime[1] = false;
    for (long long i = 2; i * i <= N; i++)
        if (isPrime[i])
            for (long long j = i * i; j <= N; j += i) isPrime[j] = false;
    vector<long long> primes;
    for (long long i = 2; i <= N; i++) if (isPrime[i]) primes.push_back(i);
    for (long long even = 4; even <= N; even += 2) {
        bool found = false;
        for (long long p : primes) {
            if (p > even / 2) break;
            if (isPrime[even - p]) { found = true; break; }
        }
        if (!found) { printf("Counterexample found: %lld\n", even); return 0; }
    }
    double s = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    printf("Goldbach verified up to %lld in %.3f s\n", N, s);
}
