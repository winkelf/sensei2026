# sensei2026

Analysis code for SENSEI Skipper-CCD data: pixel-masking studies, cluster/event
statistics ("q-statistic"), and rate plots comparing masking strategies.

## Structure

```
code/
├── cxx/        # ROOT (C++) macros only — run them from inside cxx/
│               #   helperFunctions.C (shared), multichannel.C, qWeighted.C, qWeightedSeparate.C,
│               #   LEC.C, binomialMask.C, run_thr0p1.C, rebinDistributions.C, scanCoverage.C,
│               #   q1MaskGif.C, qAppletExport.C, hduPlotter.C
├── python/     # rates.py (rates from txts/ -> results/rates/rates.pdf), qApplet.py
│               # run them from inside python/
├── configs/    # *.yaml: config.yaml (OHDUs, masks, input files), qweighted_cfg, qapplet_cfg, gif_cfg
├── scripts/    # shell launchers: run_qweighted_scan.sh, run_qweighted_separate_scan.sh, run_qapplet.sh
├── txts/       # tables written by the macros, read by python/
│   ├── lec/                      # LEC.C (rewritten 2026-09-22, same loader/mask as qWeighted)
│   ├── multichannel/r{15,20,25,30}/
│   ├── qweighted/<run>/          # qWeighted.C, one dir per (radius, w2, w3)
│   ├── qweighted_separate/<run>/
│   └── qapplet/                  # qAppletExport.C output for qApplet.py
├── root/       # .root histograms (threshold=0.1 pass), same layout as txts/
├── gifs/       # q1MaskGif.C animations
├── results/    # PDFs only
│   ├── rates/                    # rates.pdf from python/rates.py (qWeighted vs LEC)
│   ├── multichannel/r{15,..,30}/, thr0p1_rebinned/
│   ├── qweighted/<run>/, early_test_r15_w1/
│   ├── qweighted_separate/<run>/
│   ├── diagnostics/              # scanCoverage
│   └── early_exploration/        # exposure maps and q distributions from the first studies
└── logs/       # (gitignored) job logs from scripts/
```

Macros write to `../txts/<study>/`, `../results/<study>/` and `../root/<study>/`.

## Pipeline

1. `configs/config.yaml` defines which OHDUs, pixel masks (neighbor, bleed, halo,
   crosstalk, edge, hot pixel/column, LEC, etc.), and input `.root` files to use.
2. The ROOT macros in `cxx/` (`LEC.C`, `multichannel.C`, `binomialMask.C`) process
   the hit data, apply the configured masks, and write survival/rate tables to
   `txts/`, with diagnostic PDFs under `results/` and root files under `root/`.
3. `python/rates.py` reads the qWeighted and LEC tables in `txts/` and produces
   `results/rates/rates.pdf`.
4. `run_thr0p1.C`, `scanCoverage.C`, and `q1MaskGif.C` are diagnostic tools built on
   top of that pipeline: rerunning just the threshold=0.1 pass, visualizing sliding-window
   scan coverage, and animating how the q1 masking cut evolves across thresholds.

## Dependencies

- [ROOT](https://root.cern) (CERN)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
- Python: `matplotlib`, `numpy`, `scipy`
