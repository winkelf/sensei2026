# Manu

ROOT (C++) macros for a circular-window multi-electron study. For each
detector image (RUNID/LTA/OHDU) a grid of overlapping circular windows
(radius r, centres every r pixels) is laid over the sensor, and for every
window we count the reconstructed hits with exactly 1, 2, 3 or 4+ electrons
(`k1`-`k4`) and the number of live (unmasked) pixels (`N`); `q_i = k_i/N`.
The study is run for two quad selections (`good_quads`, `notbad_quads`); the
current plots all use r = 30.

## Workflow

Each dataset lives in its own folder with its own `config.yaml`; all outputs
are written inside that folder, so datasets never overwrite each other.

```
./run_dataset.sh <dataset_folder>                 # all stages
./run_dataset.sh <dataset_folder> scan lec rates  # only some stages
```

| stage | macro | output |
|---|---|---|
| `scan` | `MaskScanManye.C` | `windows` trees, standard mask (0x56FC) |
| `lec` | `MaskScanManyeLEC.C` | `windows` trees with the LEC bit (0x2000) added, in `<quads>/lec/` |
| `rates` | `ImageRates.C` | per-image live pixels and hit counts n1-n4, in `imagerates/` |
| `plots` | `PlotK2EqualsTwo.C`, `PlotQScenarios.C`, `PlotK1Independence.C` | the r = 30 PDFs |

The generation stages skip any output file that already exists, so an
interrupted run can simply be restarted. Afterwards:

```
root -l -b -q K1ClosureTest.C   # closure test of the independence method (from a dataset folder)
python3 merge_results.py        # one summary PDF per dataset and quad selection, in summary/<quads>/
```

## Files

Data generation
- `MaskScanManye.C` -- builds the per-image pixel mask from `calPixTree`,
  scans the window grid and writes one small `windows` tree per (quads,
  exposure, radius, input file): `runID, LTA, OHDU, x, y, N, k1-k4, r1-r4,
  radius, exposure, maskBits`. Windows with N = 0 are not stored.
- `MaskScanManyeLEC.C` -- the same scan with the LEC mask bit (0x2000) added,
  the reference for the q-cut scans.
- `ImageRates.C` -- per image: live pixels `Npix` and hits `n1-n4` on live
  pixels (window sums can't be used because windows overlap).

Plots (all read only the small `windows` trees, so they run in seconds)
- `PlotK2EqualsTwo.C` -- distributions of N, k and q for all windows (no cut)
  and for the selections k2 = 1, k2 = 2 and k3+ = 1 (k3+ = k3 + k4).
- `PlotQScenarios.C` -- q1-q4 distributions for no cut and k2 < 2 (log axis,
  q = 0 counted separately), and the independent per-channel q_i cut scans
  with the LEC reference lines.
- `PlotK1Independence.C` -- independence tests: the k1 (or k2) distribution of
  windows passing a tag (k2 >= 1, or k3+ = 1) against the expectation from
  independence, where each window contributes p = 1 - exp(-N lambda) (at least
  one) or N lambda exp(-N lambda) (exactly one), lambda = image rate of the
  tagging hits. Includes the ratio, chi2, and P(tag) vs the tested count.
- `K1ClosureTest.C` -- closure test of that method: in images with no 2e hit,
  fake 2e hits are placed uniformly on live pixels at the real rate and the
  same test is run (20 toys), plus a calibration page (p-values, pulls).
- `PlotWindowStats.C` -- shared helpers (`buildChain`,
  `discoverTaggedSubdirs`) used by the plotting macros.
- `PlotStyle.C` -- shared plot style (fonts, colours per charge: 1e grey,
  2e blue, 3e orange, 4e+ aqua), legends, headers and PDF page cropping.
- `PlotQ1CutScan.C` -- older q-cut scan with weighted combinations, replaced
  by `PlotQScenarios.C`; kept for reference.

Window-flag masks
- `BuildFlagMasks.C` -- writes a full copy of every input file (inputs are
  never modified) with three new mask bits in `calPixTree`'s `mask`, set on
  every pixel of a flagged r = 30 window (standard mask 0x56FC):

  | bit | flag | window flagged if |
  |---|---|---|
  | 0x100000 | A | k2 >= 2, or k2 = 1 and k1 >= 2 |
  | 0x200000 | B | k2 >= 2, or k2 = 1 and k1 >= 3 |
  | 0x400000 | C | k2 >= 2, or k2 = 1 and k1 >= 5 |

  `hitSumm` keeps its `flag` (OR of the masks of the event's pixels) and gets
  `flagA`, `flagB`, `flagC` = `flag` plus that flag's bit if any pixel of the
  event carries it. The copy also has `flagWindows` (the flagged windows) and
  `flagSummary` (per-image counts).
  `root -l -b -q 'BuildFlagMasks.C+("<raw folder>","<output folder>")'`
- `FlagSummary.C` -- how much each flag removes per dataset, quad selection
  and exposure.
- `FlagDisplay.C` -- charge map of one OHDU (white 0 e, yellow 1 e, red 5+ e)
  with the standard mask in grey and the flags (or, with `lec = true`, the
  LEC mask) in blue; writes PNG.

Animation and utilities
- `MaskAnimation.C` + `frames2gif.py` -- step-by-step frames (and a GIF) of
  how the window scan is applied to one image.
- `run_dataset.sh` -- runs the stages above for one dataset folder.
- `merge_results.py` -- merges the r = 30 PDFs into one titled, bookmarked
  PDF per dataset and quad selection (edit its `DATASETS` list for your
  folders).
- `pngs2pdf.py` -- combines PNGs into a PDF without recompression.
- `config.yaml` -- quad selections, radii and input files. **Edit
  `exposures` to point at your own copy of the data**; the paths here are
  placeholders.

## Requirements

- [ROOT](https://root.cern) with `yaml-cpp` (as the rest of this repo)
- Python 3 with Pillow (`frames2gif.py`, `pngs2pdf.py`)
- poppler-utils (`pdftocairo`, `pdfinfo`) and Ghostscript (PDF cropping and
  `merge_results.py`)

## Output layout (per dataset folder)

```
<quads>/e<exposure>/r30/      windowStatisticsManye_<quads>_<n>.root and the plot PDFs
<quads>/lec/e<exposure>/r30/  LEC-mask windows
imagerates/e<exposure>/       imageRates_<n>.root
```

`<n>` is the position of the input file in the config's list for that
exposure; it is also used to match windows to image rates, because run IDs
repeat between data campaigns.

## Config format

```yaml
good_quads:
  <LTA>: [<OHDU>, ...]   # quads considered good, per LTA
notbad_quads:
  <LTA>: [<OHDU>, ...]   # quads considered usable (superset of good_quads)
radii: [30]              # window radii (pixels) to scan
exposures:
  <exposure>:            # e.g. 216000
    - /path/to/file1.root
    - ...
```

LTAs/OHDUs not listed under a given selection are excluded from it.
