"""n-electron rates from the txt tables: qWeighted.C (single combined
statistic q = (k1 + w2*k2 + w3*(k3+k4+k5))/N, one shared mask) vs LEC.

Both come from the same image loader (loadImage3: events at their barycenter
pixel, e truncated), the same event definition (1e / 2e / 3+4+5e), the same
GOOD_QUADS set and the same standard SENSEI base mask (0x56FC). The only
difference is the extra mask on top: the sliding-window combined-q mask for
qWeighted, the static LEC bit (8192) for LEC.

Page 1: per-channel rates for w2 = w3 = 1 (r=15) vs LEC and the base-mask-only
        rate, with 1-sigma Poisson bands sqrt(n)/N_pix.
Pages 2+: one page per (w2, w3) weight pair at r=15 (threshold scan
        0.0012-0.02): surviving events, unmasked pixels and rate vs the q cut.
        LEC applies one fixed static mask (no sliding-window cut), so its
        numbers are threshold-independent and drawn as flat lines.

Inputs:  ../txts/qweighted/r15_w2-<w2>_w3-<w3>/, ../txts/lec/
Output:  ../results/rates/rates.pdf
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

radius = 15
pairs = [(1, 1), (10, 10), (10, 100), (100, 10), (100, 1000), (1000, 100),
         (0.1, 1), (0.1, 0.1), (0.01, 1), (0.01, 0.1)]

QW_DIR  = "../txts/qweighted"
LEC_DIR = "../txts/lec"
pdfName = "../results/rates/rates.pdf"

C_QW   = "#2a78d6"
C_LEC  = "#eb6834"
C_BASE = "#8a8a85"
INK    = "#1f1f1e"
MUTED  = "#6b6b66"
GRID   = "#e4e3dc"

LIGHT_RC = {
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
}

lecSurvElec = np.atleast_2d(np.loadtxt(f"{LEC_DIR}/LEC_survElec.txt"))
lecSurvPix  = np.atleast_2d(np.loadtxt(f"{LEC_DIR}/LEC_survPix.txt"))
lecRates    = np.atleast_2d(np.loadtxt(f"{LEC_DIR}/LEC_rates.txt"))


def hband(ax, n, pix, color, label, ls):
    r, e = n/pix, np.sqrt(n)/pix
    ax.axhspan(r - e, r + e, color=color, alpha=0.15, lw=0)
    ax.axhline(r, color=color, lw=2, ls=ls, label=label)


def comparison_page(pdf):
    """w2 = w3 = 1 per-channel rates vs LEC and base-mask-only."""
    d = f"{QW_DIR}/r{radius}_w2-1_w3-1"
    survElec = np.loadtxt(f"{d}/qweighted_survElec.txt")
    survPix  = np.loadtxt(f"{d}/qweighted_survPix.txt")

    # Float-accumulated threshold loop yields 0.002 twice; keep the first.
    _, keep = np.unique(np.round(survElec[:, 0], 6), return_index=True)
    survElec, survPix = survElec[keep], survPix[keep]

    thr  = survElec[:, 0]
    nPix = survPix[:, 1]

    lec = np.loadtxt(f"{LEC_DIR}/LEC_perImage.txt")
    # columns: lta ohdu runId n1_base n2_base n345_base nPix_base n1_LEC n2_LEC n345_LEC nPix_LEC
    base_n,  base_pix = lec[:, 3:6].sum(axis=0), lec[:, 6].sum()
    lec_n,   lec_pix  = lec[:, 7:10].sum(axis=0), lec[:, 10].sum()

    channels = ["1 electron", "2 electrons", "3+4+5 electrons"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"n-electron rates per channel: combined-q mask ($w_2=w_3=1$, r={radius}) vs LEC\n"
                 "same loader, event definition, GOOD_QUADS and base mask (0x56FC)",
                 color=INK, fontsize=12)

    for c, ax in enumerate(axes.flat[:3]):
        n = survElec[:, c+1]
        r, e = n/nPix, np.sqrt(n)/nPix

        hband(ax, base_n[c], base_pix, C_BASE, "Base mask only", ":")
        hband(ax, lec_n[c],  lec_pix,  C_LEC,  "Base + LEC", "--")
        ax.fill_between(thr, r - e, r + e, color=C_QW, alpha=0.2, lw=0)
        ax.plot(thr, r, "-o", color=C_QW, lw=2, ms=4, label="Base + combined-q mask")

        ax.set_xscale("log")
        ax.set_title(f"{channels[c]} rate", loc="left")
        ax.set_xlabel("Cut in $q=(k_1+k_2+k_3+k_4+k_5)/N$")
        ax.set_ylabel("Rate [events / unmasked pixel]")
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.text(0.98, 0.05,
                f"LEC: {int(lec_n[c])} ev   base: {int(base_n[c])} ev   "
                f"q-mask @0.02: {int(n[-1])} ev",
                transform=ax.transAxes, ha="right", color=MUTED, fontsize=8)

    ax = axes.flat[3]
    ax.axhline(base_pix, color=C_BASE, lw=2, ls=":", label="Base mask only")
    ax.axhline(lec_pix,  color=C_LEC,  lw=2, ls="--", label="Base + LEC")
    ax.plot(thr, nPix, "-o", color=C_QW, lw=2, ms=4, label="Base + combined-q mask")
    ax.set_xscale("log")
    ax.set_title("Unmasked pixels (exposure)", loc="left")
    ax.set_xlabel("Cut in $q=(k_1+k_2+k_3+k_4+k_5)/N$")
    ax.set_ylabel("Unmasked pixels")
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)

    ticks = [0.0015, 0.002, 0.003, 0.005, 0.01, 0.02]
    for ax in axes.flat:
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.minorticks_off()

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.925), labelcolor=INK)

    fig.tight_layout(rect=(0, 0, 1, 0.89))
    pdf.savefig(fig)
    plt.close(fig)


def weight_pair_page(pdf, w2, w3):
    """Surviving events, unmasked pixels and rate vs q cut for one (w2, w3)."""
    d = f"{QW_DIR}/r{radius}_w2-{w2}_w3-{w3}"
    survElec = np.loadtxt(f"{d}/qweighted_survElec.txt")
    survPix  = np.loadtxt(f"{d}/qweighted_survPix.txt")
    rates    = np.loadtxt(f"{d}/qweighted_rates.txt")

    x = survElec[:, 0]

    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax0 = fig.add_subplot(gs[0, 0])  # Top left
    ax1 = fig.add_subplot(gs[1, 0])  # Bottom left
    ax2 = fig.add_subplot(gs[:, 1])  # Right (spans both rows)

    # Plot 1: surviving events, one curve per channel.
    y1, y2, y345 = survElec[:, 1], survElec[:, 2], survElec[:, 3]
    y1LEC   = np.full_like(x, lecSurvElec[0, 1])
    y2LEC   = np.full_like(x, lecSurvElec[0, 2])
    y345LEC = np.full_like(x, lecSurvElec[0, 3])

    ax0.plot(x, y1,   "-o", label="1 electron count", color="blue", markersize=3)
    ax0.plot(x, y2,   "-o", label="2 electrons count", color="green", markersize=3)
    ax0.plot(x, y345, "-o", label="3+4+5 electrons count", color="orange", markersize=3)
    ax0.plot(x, y1LEC,   "--", label="1 electron count (LEC)", color="blue")
    ax0.plot(x, y2LEC,   "--", label="2 electrons count (LEC)", color="green")
    ax0.plot(x, y345LEC, "--", label="3+4+5 electrons count (LEC)", color="orange")
    ax0.set_title("Number of surviving n-electron events vs q cut")
    ax0.set_xlabel('Cut in $q=(k_1+w_2k_2+w_3(k_3+k_4+k_5))/N$')
    ax0.set_ylabel("Number of surviving n-electron events")
    ax0.grid(True, linestyle=':')
    ax0.legend(fontsize=8)

    # Plot 2: unmasked pixels (single shared mask).
    y_pix = survPix[:, 1]
    yLEC_pix = np.full_like(x, lecSurvPix[0, 1])

    ax1.plot(x, y_pix, "-o", label="Unmasked pixels", color="blue", markersize=3)
    ax1.plot(x, yLEC_pix, "--", label="Unmasked pixels (LEC)", color="yellow")
    ax1.set_title("Number of unmasked pixels vs q cut")
    ax1.set_xlabel('Cut in $q$')
    ax1.set_ylabel("Number of unmasked pixels")
    ax1.grid(True, linestyle=':')
    ax1.legend(fontsize=8)

    # Plot 3: rates, one curve per channel.
    rate1, rate2, rate345 = rates[:, 1], rates[:, 2], rates[:, 3]
    rate1LEC   = np.full_like(x, lecRates[0, 1])
    rate2LEC   = np.full_like(x, lecRates[0, 2])
    rate345LEC = np.full_like(x, lecRates[0, 3])

    ax2.plot(x, rate1,   "-o", label="1 electron rate", color="blue", markersize=3)
    ax2.plot(x, rate2,   "-o", label="2 electrons rate", color="green", markersize=3)
    ax2.plot(x, rate345, "-o", label="3+4+5 electrons rate", color="orange", markersize=3)
    ax2.plot(x, rate1LEC,   "--", label="1 electron rate (LEC)", color="blue")
    ax2.plot(x, rate2LEC,   "--", label="2 electrons rate (LEC)", color="green")
    ax2.plot(x, rate345LEC, "--", label="3+4+5 electrons rate (LEC)", color="orange")
    ax2.set_title("Rates vs q cut")
    ax2.set_xlabel('Cut in $q$')
    ax2.set_ylabel("n-electron rate")
    ax2.set_yscale('symlog', linthresh=1e-7)
    ax2.set_ylim(bottom=0)
    ax2.grid(True, linestyle=':')
    ax2.legend(fontsize=8)

    fig.suptitle(f"qWeighted (combined statistic, shared mask), r={radius}, w2={w2}, w3={w3}")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages(pdfName) as pdf:
    with plt.rc_context(LIGHT_RC):
        comparison_page(pdf)
    with plt.style.context('dark_background'):
        for w2, w3 in pairs:
            weight_pair_page(pdf, w2, w3)

print(f"Generated {pdfName}")
