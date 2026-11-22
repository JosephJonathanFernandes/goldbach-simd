# goldbach-simd

High-throughput C++ Goldbach conjecture verifier using a cache-blocked segmented bucket sieve, AVX2 bit-parallel checks, multithreading, and checkpointing.

## Highlights

- **Cache-blocked segmented sieve** with $O(\sqrt{N})$ memory (~12 MB total across 12 threads at $10^{12}$).
- **Pre-sieved periodic bit patterns** eliminating all odd primes up to 43 before crossing off composites.
- **L1-threshold bucket sieve** ($p \ge 2^{18}$) + **contiguous chunk batching** to amortize 64-bit divisions.
- **4-prime fused AVX2 SIMD kernel** checking 256 even numbers per iteration in YMM registers (with SSE2 and scalar fallbacks).
- **Deterministic 64-bit Miller-Rabin fallback** so no even number can ever be falsely verified.
- **Built-in $\pi(N)$ self-check** against known prime counts up to $10^{12}$ and **checkpoint/resume** support.

## Benchmarks (12 threads)

| $N$ | Baseline (`baseline.cpp`) | `goldbach.cpp` (12T) | Throughput | Largest min $p$ ($n = p + q$) | $\pi(N)$ |
|---:|---:|---:|---:|---|---:|
| $10^8$ | 2.315 s | **0.093 s** | 539 M evens/s | `1093` ($60{,}119{,}912$) | `5,761,455` |
| $10^9$ | 27.261 s | **~0.16 s** (0.908 s 1T) | — | `1789` ($721{,}013{,}438$) | `50,847,534` |
| $10^{10}$ | *(1.25 GB RAM)* | **1.154 s** | **4.33 B evens/s** | `1877` ($7{,}473{,}202{,}036$) | `455,052,511` |
| $10^{11}$ | *(12.5 GB RAM)* | **13.827 s** | **3.62 B evens/s** | `2803` ($35{,}884{,}080{,}836$) | `4,118,054,813` |
| $10^{12}$ | *(125 GB RAM)* | **200.813 s** | **2.49 B evens/s** | `3457` ($599{,}533{,}546{,}358$) | `37,607,912,018` |

## Build & Run

```bash
# Linux / macOS / MSYS2
g++ -O3 -march=native -std=c++17 -pthread goldbach.cpp -o goldbach

# Windows (MinGW)
g++ -O3 -march=native -std=c++14 goldbach.cpp -o goldbach.exe
```

Usage:

```bash
./goldbach [N] [threads] [checkpoint_file]

# Examples
./goldbach 1e10
./goldbach 1e12 12 goldbach_1e12.ckpt
```
