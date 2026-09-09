# Manu

ROOT (C++) macros for a circular-window multi-electron rate study: for each
detector image (RUNID/LTA/OHDU), a grid of overlapping circular windows is
laid over the sensor, and for every window we count how many reconstructed
hits have exactly 1, 2, 3, or 4+ electrons (`k1`-`k4`), normalized by the
number of valid (unmasked) background pixels in that window (`N`). The study
is run separately for two quad-quality selections (`good_quads`,
`notbad_quads`), and across several window radii, so you can compare which
mask size gives the most stable rate estimate.

## Files

- `MaskScanManye.C` -- the data-generation macro. Reads `hitSumm`/`calPixTree`
  from each input file, builds the per-detector pixel mask, scans the window
  grid, and writes one small `windows` TTree per (category, radius, input
  file) with branches `runID, LTA, OHDU, x, y, N, k1-k4, r1-r4, radius`.
- `PlotWindowStats.C` -- the plotting macro. Chains together the `windows`
  trees MaskScanManye.C wrote and produces summary PDFs. It only reads the
  small per-window trees, not the original hit data, so it's fast and can be
  re-run any time without redoing the analysis.
- `config.yaml` -- defines the quad selections, radii to scan, and input
  files. **Edit `input_files` to point at your own local copy of the data**
  before running -- the paths shipped here are placeholders.

## Requirements

- [ROOT](https://root.cern) (with `yaml-cpp` available -- same as the rest of
  this repo, see the top-level README)

## How to run

From inside this folder:

```
root -l -b -q MaskScanManye.C
root -l -b -q PlotWindowStats.C
```

`MaskScanManye.C` will take a while for a full radius sweep over many files
(it's O(categories x radii x files)); `PlotWindowStats.C` runs in seconds
since it only re-reads the small output trees.

## Output layout

Everything is written under `<category>/r<radius>/`, e.g.:

```
goodquads/
  r15/  windowStatisticsManye_goodquads_1.root ... _60.root
        k_statistics_summary_goodquads_r15.pdf
        q_statistics_summary_goodquads_r15.pdf
  r20/  ...
  r25/  ...
  r30/  ...
notbadquads/
  r15/  ...
  ...
```

- `windowStatisticsManye_<category>_<n>.root` -- one per input file, contains
  the `windows` TTree (raw per-window data) and quick-look ratio histograms.
- `k_statistics_summary_<category>_r<radius>.pdf` -- k1-k4 (raw counts) and N
  distributions across all windows/files for that category and radius, plus
  a stats box with totals.
- `q_statistics_summary_<category>_r<radius>.pdf` -- same, but for
  q1-q4 = k1-k4 / N (per-window rates), plus N.

To decide which radius is optimal, compare the summary PDFs across
`r15/r20/r25/r30` for each category -- look at how the q distributions and
their statistical spread (fewer windows but larger N vs. more windows but
noisier per-window counts) change with radius.

## Config format

```yaml
good_quads:
  <LTA>: [<OHDU>, ...]   # quads considered good, per LTA
notbad_quads:
  <LTA>: [<OHDU>, ...]   # quads considered usable (superset of good_quads)
radii: [15, 20, 25, 30]  # window radii (pixels) to scan
input_files:
  - /path/to/file1.root
  - ...
```

LTAs/OHDUs not listed under a given category are excluded from that
category's analysis entirely.
