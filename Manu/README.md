# Manu

ROOT (C++) macros for a circular-window multi-electron rate study: for each
detector image (RUNID/LTA/OHDU), a grid of overlapping circular windows is
laid over the sensor, and for every window we count how many reconstructed
hits have exactly 1, 2, 3, or 4+ electrons (`k1`-`k4`), normalized by the
number of valid (unmasked) background pixels in that window (`N`). The study
is run separately for two quad-quality selections (`good_quads`,
`notbad_quads`), across several window radii (to compare which mask size
gives the most stable rate estimate), and across exposures (to see whether
mask performance changes with exposure time).

## Files

- `MaskScanManye.C` -- the data-generation macro. Reads `hitSumm`/`calPixTree`
  from each input file, builds the per-detector pixel mask, scans the window
  grid, and writes one small `windows` TTree per (category, exposure, radius,
  input file) with branches `runID, LTA, OHDU, x, y, N, k1-k4, r1-r4, radius,
  exposure`.
- `PlotWindowStats.C` -- the plotting macro. Chains together the `windows`
  trees MaskScanManye.C wrote and produces summary PDFs. It only reads the
  small per-window trees, not the original hit data, so it's fast and can be
  re-run any time without redoing the analysis. It auto-discovers every
  `e<exposure>/r<radius>` combination present on disk, so it doesn't need to
  read `config.yaml` at all.
- `MaskScanManyeLEC.C` -- runs the same category x exposure x radius x file
  scan as `MaskScanManye.C`, but with the LEC ("Local Event Cut") pixel-mask
  bit (0x2000) added on top of the standard mask, writing to
  `<category>/lec/e<exposure>/r<radius>/`. This gives a single (k1-k4, N)
  reference per combination to compare against a q1 cut on the standard-mask
  output -- see `PlotQ1CutScan.C`. `#include`s `MaskScanManye.C` to reuse its
  types/helpers/`windowStudy()`, and only defines its own `MaskScanManyeLEC()`
  entry point, so running it doesn't touch or require the standard-mask
  output.
- `PlotQ1CutScan.C` -- for each (category, exposure, radius), scans a q1 =
  k1/N cut from 0 up to `min(0.02, highest observed q1)`, keeping only
  windows with q1 < cut at each step, and plots k1-k4, N, and k_i/N vs. the
  cut (3 PDFs per combination). Each plot overlays the single LEC reference
  value (dashed line) from `MaskScanManyeLEC.C`'s output, so you can see
  directly whether tightening the q1 cut can match what the LEC mask
  achieves. `#include`s `PlotWindowStats.C` to reuse its helpers.
- `config.yaml` -- defines the quad selections, radii to scan, and input
  files grouped by exposure. **Edit `exposures` to point at your own local
  copy of the data** before running -- the paths shipped here are
  placeholders.

## Requirements

- [ROOT](https://root.cern) (with `yaml-cpp` available -- same as the rest of
  this repo, see the top-level README)

## How to run

From inside this folder:

```
root -l -b -q MaskScanManye.C
root -l -b -q PlotWindowStats.C
```

and, if you also want the LEC-mask comparison plots:

```
root -l -b -q MaskScanManyeLEC.C
root -l -b -q PlotQ1CutScan.C
```

`MaskScanManye.C` and `MaskScanManyeLEC.C` will each take a while for a full
sweep over many files (it's O(categories x exposures x radii x files));
`PlotWindowStats.C` and `PlotQ1CutScan.C` run in seconds since they only
re-read the small output trees. If a run is interrupted partway through,
just re-run the same command -- both generation macros skip any output file
that already exists, so nothing gets recomputed.

## Output layout

Everything is written under `<category>/e<exposure>/r<radius>/`, e.g.:

```
goodquads/
  e72000/
    r15/  windowStatisticsManye_goodquads_1.root ... _60.root
          k_statistics_summary_goodquads_e72000_r15.pdf
          q_statistics_summary_goodquads_e72000_r15.pdf
          q1cut_k_scan_goodquads_e72000_r15.pdf
          q1cut_N_scan_goodquads_e72000_r15.pdf
          q1cut_ratio_scan_goodquads_e72000_r15.pdf
    r20/  ...
    r25/  ...
    r30/  ...
  e216000/
    r15/  ...
    ...
  lec/
    e72000/
      r15/  windowStatisticsManye_goodquads_1.root ... _60.root  (LEC mask)
      r20/  ...
      ...
    e216000/
      ...
notbadquads/
  e72000/
    ...
  e216000/
    ...
  lec/
    ...
```

- `windowStatisticsManye_<category>_<n>.root` -- one per input file, contains
  the `windows` TTree (raw per-window data, including its own `exposure`,
  `radius` and `maskBits` branches) and quick-look ratio histograms. Under
  `lec/`, these were produced with the LEC mask bit added instead of the
  standard mask.
- `k_statistics_summary_<category>_e<exposure>_r<radius>.pdf` -- k1-k4 (raw
  counts) and N distributions across all windows/files for that category,
  exposure and radius, plus a stats box with totals.
- `q_statistics_summary_<category>_e<exposure>_r<radius>.pdf` -- same, but
  for q1-q4 = k1-k4 / N (per-window rates), plus N.
- `q1cut_k_scan_<category>_e<exposure>_r<radius>.pdf`,
  `q1cut_N_scan_..._.pdf`, `q1cut_ratio_scan_..._.pdf` -- k1-k4, N, and
  k_i/N as a function of a q1 = k1/N cut (windows with q1 above the cut are
  dropped), each with the corresponding LEC-mask reference value overlaid as
  a dashed line. Produced by `PlotQ1CutScan.C` from the standard-mask output
  plus the single reference computed once from `lec/`.

To decide which radius is optimal, compare the summary PDFs across
`r15/r20/r25/r30` for a fixed category and exposure -- look at how the q
distributions and their statistical spread (fewer windows but larger N vs.
more windows but noisier per-window counts) change with radius. To check
whether mask performance changes with exposure, compare across
`e72000/e216000` at a fixed radius instead.

## Config format

```yaml
good_quads:
  <LTA>: [<OHDU>, ...]   # quads considered good, per LTA
notbad_quads:
  <LTA>: [<OHDU>, ...]   # quads considered usable (superset of good_quads)
radii: [15, 20, 25, 30]  # window radii (pixels) to scan
exposures:
  <exposure>:            # e.g. 72000, 216000
    - /path/to/file1.root
    - ...
```

LTAs/OHDUs not listed under a given category are excluded from that
category's analysis entirely.
