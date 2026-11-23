"""Generate publication-grade figures for goldbach-simd RESEARCH.md and README.md."""
import os
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("figures", exist_ok=True)

plt.rcParams.update({
    "figure.facecolor": "#0d1117",
    "axes.facecolor": "#161b22",
    "axes.edgecolor": "#30363d",
    "axes.labelcolor": "#e6edf3",
    "text.color": "#e6edf3",
    "xtick.color": "#8b949e",
    "ytick.color": "#8b949e",
    "grid.color": "#21262d",
    "grid.linestyle": "--",
    "grid.alpha": 0.8,
    "font.size": 11,
    "axes.titleweight": "bold",
    "axes.titlesize": 13,
})

# ---------------------------------------------------------------------------
# Figure 1: Exact Global Carry-Optimal Goldbach Phase Transition (1e8 vs 1e9)
# ---------------------------------------------------------------------------
w_1e8 = np.arange(1, 12)
fail_1e8 = np.array([
    100.0, 96.91358, 89.59627, 80.91863, 65.24406,
    46.52921, 25.86645, 9.59283, 1.85050, 0.09584, 0.00034
])

w_1e9 = np.arange(1, 13)
fail_1e9 = np.array([
    100.0, 97.53695, 91.46141, 84.20632, 70.90963,
    53.94064, 33.74183, 15.07480, 4.02328, 0.38307, 0.00605, 100.0 / 50777720.0
])

fig, ax = plt.subplots(figsize=(9, 5.2), dpi=160)
ax.semilogy(w_1e8, fail_1e8, "o-", color="#58a6ff", lw=2.3, ms=7,
            label=r"$N = 10^8$ (50M evens — 0 exceptions for $w \geq 12$)")
ax.semilogy(w_1e9, fail_1e9, "s-", color="#f78166", lw=2.3, ms=7,
            label=r"$N = 10^9$ (500M evens — 0 exceptions for $w \geq 13$)")

ax.annotate(
    "Unique w=12 needle\nn = 743,309,774 (0x2C4E01CE)\n1 / 50,777,720",
    xy=(12, 100.0 / 50777720.0), xytext=(8.6, 4e-5),
    arrowprops=dict(arrowstyle="->", color="#ffa657", lw=1.6),
    fontsize=9.5, color="#ffa657",
    bbox=dict(boxstyle="round,pad=0.35", fc="#1f242c", ec="#ffa657", lw=1)
)
ax.axvspan(12.5, 15.2, color="#238636", alpha=0.14)
ax.text(13.85, 0.12, "Zero Exceptions\n(378.8M evens\nfor w = 13..28)",
        ha="center", va="center", color="#3fb950", fontsize=10, fontweight="bold")

ax.set_xlim(0.6, 15.2)
ax.set_ylim(5e-7, 250)
ax.set_xticks(range(1, 16))
ax.set_xlabel(r"Binary Weight $w = \mathrm{popcount}(n - 2)$")
ax.set_ylabel("Global Exception Rate (%) — Log Scale")
ax.set_title(r"Exact Global Carry-Optimal Goldbach ($\Delta_c^*(n) > 0$) Phase Transition")
ax.grid(True, which="both")
ax.legend(loc="lower left", facecolor="#0d1117", edgecolor="#30363d")
plt.tight_layout()
plt.savefig("figures/carry_phase_transition.png")
plt.close()

# ---------------------------------------------------------------------------
# Figure 2: Carry Excess Tail Decay (Small-Witness 2^-k vs Exact Global)
# ---------------------------------------------------------------------------
dc_small = np.arange(0, 16)
cnt_small = np.array([
    398669118, 64054308, 24402075, 8517747, 2964012, 919796,
    292541, 102108, 40524, 19077, 9328, 4666, 2359, 1150, 597, 286
], dtype=float)
pct_small = 100.0 * cnt_small / 499999998.0

dc_glob = np.arange(0, 10)
cnt_glob = np.array([
    49489773, 438377, 59418, 10088, 1918, 321, 73, 21, 8, 1
], dtype=float)
pct_glob = 100.0 * cnt_glob / 49999998.0

# Reference 2^(-Delta_c) line anchored at Delta_c = 8
ref_x = np.arange(7, 16)
ref_y = pct_small[8] * (2.0 ** (8 - ref_x))

fig, ax = plt.subplots(figsize=(9, 5.2), dpi=160)
ax.semilogy(dc_small, pct_small, "o-", color="#79c0ff", lw=2.3, ms=6.5,
            label=r"Small-witness $\Delta_c(n)$ ($p < 16384$, $N = 10^9$)")
ax.semilogy(ref_x, ref_y, "--", color="#ffa657", lw=2.0,
            label=r"Exact geometric law $\propto 2^{-\Delta_c}$ (zero-run borrow cascade)")
ax.semilogy(dc_glob, pct_glob, "D-", color="#d2a8ff", lw=2.3, ms=6.5,
            label=r"Exact global $\Delta_c^*(n)$ (all $p \leq n/2$, $N = 10^8$)")

ax.annotate(
    "Unique global max: Δc* = 9\nn = 8,388,610 = 2²³ + 2",
    xy=(9, pct_glob[9]), xytext=(9.6, 5e-6),
    arrowprops=dict(arrowstyle="->", color="#d2a8ff", lw=1.6),
    fontsize=9.5, color="#d2a8ff",
    bbox=dict(boxstyle="round,pad=0.35", fc="#1f242c", ec="#d2a8ff", lw=1)
)

ax.set_xlim(-0.4, 15.5)
ax.set_xlabel(r"Carry Excess $\Delta_c(n) = \mathrm{carries}(p, q) - v_2(n)$")
ax.set_ylabel("Share of Even Integers (%) — Log Scale")
ax.set_title("Goldbach Binary Carry Excess Distribution: Small-Witness vs. Exact Global")
ax.grid(True, which="both")
ax.legend(loc="upper right", facecolor="#0d1117", edgecolor="#30363d")
plt.tight_layout()
plt.savefig("figures/carry_tail_decay.png")
plt.close()

# ---------------------------------------------------------------------------
# Figure 3: Minimal Goldbach Witness Frequencies & Chebyshev Mod-6 Inversion
# ---------------------------------------------------------------------------
primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
freq = [9.1004, 8.5522, 8.0041, 7.0165, 6.4178, 5.5387, 5.3590, 4.5206, 3.9865, 4.0373]
colors = [
    "#8b949e" if p == 3 else ("#f78166" if p % 6 == 5 else "#58a6ff")
    for p in primes
]

fig, ax = plt.subplots(figsize=(9, 5.0), dpi=160)
bars = ax.bar([str(p) for p in primes], freq, color=colors, edgecolor="#30363d", width=0.66)
for bar, f in zip(bars, freq):
    ax.text(bar.get_x() + bar.get_width() / 2, f + 0.14, f"{f:.2f}%",
            ha="center", va="bottom", fontsize=9, color="#e6edf3")

ax.annotate(
    "Chebyshev mod-6 inversion!\np=31 (4.04%) > p=29 (3.99%)",
    xy=(9, 4.04), xytext=(6.0, 6.9),
    arrowprops=dict(arrowstyle="->", color="#3fb950", lw=1.8),
    fontsize=9.5, color="#3fb950", fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.35", fc="#1f242c", ec="#3fb950", lw=1)
)

# Legend patches for mod-6 classes
from matplotlib.patches import Patch
legend_elems = [
    Patch(facecolor="#8b949e", label=r"$p = 3$"),
    Patch(facecolor="#f78166", label=r"$p \equiv 5\ (\mathrm{mod}\ 6)$"),
    Patch(facecolor="#58a6ff", label=r"$p \equiv 1\ (\mathrm{mod}\ 6)$"),
]
ax.legend(handles=legend_elems, loc="upper right", facecolor="#0d1117", edgecolor="#30363d")

ax.text(0.02, 0.94,
        r"$N = 10^{10}$ (5B evens)  |  $H(p_{\min}) = 4.8616$ bits  |  $I(X_t; X_{t+1}) = 1.5872 > \log_2 3$ bits",
        transform=ax.transAxes, fontsize=9.2, color="#c9d1d9", va="top")

ax.set_ylim(0, 10.6)
ax.set_xlabel(r"Minimal Goldbach Witness Prime $p_{\min}(n)$")
ax.set_ylabel("Frequency across 5 Billion Evens (%)")
ax.set_title(r"Minimal Goldbach Witness Distribution ($N = 10^{10}$) & Mod-6 Inversion")
ax.grid(True, axis="y")
plt.tight_layout()
plt.savefig("figures/witness_entropy_inversion.png")
plt.close()

# ---------------------------------------------------------------------------
# Figure 4: Engineering Throughput Progression (Baseline -> goldbach-simd)
# ---------------------------------------------------------------------------
stages = [
    "Baseline\nvector<bool>\n(1T)",
    "R1: Segmented\n+ 64-bit Check\n(12T)",
    "R2: Pre-sieve\n+ L1 Block\n(12T)",
    "R3: Batch Div\n+ Checkpoints\n(12T)",
    "R4–5: Fused\n4-Prime AVX2\n(12T)",
    "R6–7: L1-Bucket\n+ AVX2 Pre-sieve\n(12T)",
]
throughput = [18.3, 748.0, 2765.0, 3397.0, 3702.0, 4332.8]
bar_colors = ["#6e7681", "#388bfd", "#58a6ff", "#79c0ff", "#a5d6ff", "#3fb950"]

fig, ax = plt.subplots(figsize=(9.2, 5.0), dpi=160)
bars = ax.bar(stages, throughput, color=bar_colors, edgecolor="#30363d", width=0.6)
for bar, val in zip(bars, throughput):
    lbl = f"{val/1000:.2f}B/s" if val >= 1000 else f"{val:.0f}M/s"
    ax.text(bar.get_x() + bar.get_width() / 2, val + 95, lbl,
            ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="#e6edf3")

ax.annotate(
    "237× faster than baseline!\n10¹² verified in 200.8 s (~12 MB RAM)",
    xy=(5, 4332.8), xytext=(2.2, 4150),
    arrowprops=dict(arrowstyle="->", color="#3fb950", lw=1.8),
    fontsize=9.5, color="#3fb950", fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.35", fc="#1f242c", ec="#3fb950", lw=1)
)

ax.set_ylim(0, 5000)
ax.set_ylabel("Peak Verification Throughput (Million Evens / sec)")
ax.set_title("goldbach-simd Optimization Progression (18M → 4.33B Evens/sec)")
ax.grid(True, axis="y")
plt.tight_layout()
plt.savefig("figures/benchmark_scaling.png")
plt.close()

print("Generated all 4 figures in figures/")
