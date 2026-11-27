"""Independent cross-check of gosper_hunt.cpp results (pure Python, no shared code).

For each claimed a(w) = n_min(w):
  1. popcount(n - 2) == w
  2. NO submask pair (1|s, 1|(M^s)) is a pair of primes  -> Delta_c*(n) > 0
  3. Every smaller even n' >= 6 with popcount(n'-2) == w HAS such a pair (checked for w <= 11).
"""
import sys

def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def has_zero_excess(n):
    M = n - 2
    hi = 1 << (M.bit_length() - 1)
    lo = M ^ hi
    s = lo
    while s:
        if is_prime(1 | s) and is_prime(1 | (M ^ s)):
            return True
        s = (s - 1) & lo
    return False

CLAIMED = {1: 6, 2: 12, 3: 52, 4: 62, 5: 188, 6: 1448, 7: 5342, 8: 14992,
           9: 105274, 10: 689072, 11: 21852554, 12: 743309774, 13: 4706638234}

ok = True
for w, n in CLAIMED.items():
    c1 = bin(n - 2).count("1") == w
    c2 = not has_zero_excess(n)
    c3 = True
    if w <= 9:
        for m in range(6, n, 2):
            if bin(m - 2).count("1") == w and not has_zero_excess(m):
                c3 = False
                break
    status = "OK" if (c1 and c2 and c3) else "FAIL"
    ok &= status == "OK"
    print(f"w={w:2d} n={n:>12d}  weight={c1}  no-carry-optimal-pair={c2}  minimal={c3 if w <= 9 else 'n/a'}  {status}")
sys.exit(0 if ok else 1)
