"""Short slide deck summarizing the window-radius scan (multichannel.C) and
the LEC masking comparison. Assumes the plotting scripts (rates.py,
rates_by_radius.py) and the ROOT macros (multichannel.C, LEC.C) have already
been run, and that their outputs have been rasterized to PNG under
/tmp/deck_assets/ (see the accompanying pdftoppm calls).
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

plt.style.use('dark_background')
FIGSIZE = (13.33, 7.5)  # 16:9 slide

assets = "/tmp/deck_assets"


def text_slide(pdf, title, lines, fontsize=20):
    fig = plt.figure(figsize=FIGSIZE)
    fig.text(0.08, 0.85, title, fontsize=32, fontweight="bold", va="top")
    y = 0.65
    for line in lines:
        fig.text(0.1, y, line, fontsize=fontsize, va="top", wrap=True)
        y -= 0.11
    pdf.savefig(fig)
    plt.close(fig)


def image_slide(pdf, title, image_paths, captions):
    fig = plt.figure(figsize=FIGSIZE)
    fig.text(0.5, 0.95, title, fontsize=24, fontweight="bold", ha="center", va="top")
    n = len(image_paths)
    for i, (path, caption) in enumerate(zip(image_paths, captions)):
        ax = fig.add_axes([0.03 + i * (0.94 / n), 0.08, 0.94 / n - 0.02, 0.78])
        ax.imshow(mpimg.imread(path))
        ax.axis("off")
        ax.set_title(caption, fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages("../results/presentation.pdf") as pdf:

    text_slide(pdf, "SENSEI: Multi-Electron Rate\nWindow-Radius Scan", [
        "Sliding-window multi-electron rate study (multichannel.C)",
        "Comparing window radii r = 15, 20, 25, 30 px",
        "Reference: Local Event Cut (LEC) global masking rate",
    ])

    text_slide(pdf, "Method", [
        "For each radius, a grid of overlapping circular windows is scanned",
        "over the sensor; k1-k4 (n-electron hit counts) are normalized by",
        "the number of unmasked background pixels N in that window.",
        "",
        "Rate vs. q-cut (k/N) is computed for 30 threshold points",
        "(0 to 0.01, step 1/3000) at each radius.",
        "",
        "LEC applies a fixed global mask (not window-based), giving a",
        "single reference rate independent of the q-cut.",
    ])

    image_slide(
        pdf,
        "Rate vs. cut, across window radii",
        [f"{assets}/rates_by_radius-1.png"],
        ["r=15/20/25/30 vs. LEC reference"],
    )

    image_slide(
        pdf,
        "Masking strategy comparison at r=30",
        [f"{assets}/rates-1.png"],
        ["With vs. without LEC masking"],
    )

    image_slide(
        pdf,
        "q-statistic (1-electron) distribution: r=30 vs. r=15",
        [f"{assets}/q1e_r30-1_crop.png", f"{assets}/q1e_r15-1_crop.png"],
        ["r=30 (stable)", "r=15 (outlier)"],
    )

    text_slide(pdf, "Findings", [
        "r=20, 25 and 30 agree closely: rate onset near cut~0.0007-0.001,",
        "plateauing at a consistent rate for all n-electron categories.",
        "",
        "r=15 is a clear outlier: onset shifted right by ~0.0007-0.001,",
        "and plateaus at a lower rate across 1e/2e/3e/4e -- consistent",
        "with a window too small to give a stable per-window estimate.",
        "",
        "Recommendation: use r=20-30 for the stable rate estimate;",
        "treat r=15 as below the useful window-size range.",
    ])

print("Generated ../results/presentation.pdf")
