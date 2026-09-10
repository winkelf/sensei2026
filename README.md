# sensei2026

Analysis code for SENSEI Skipper-CCD data: pixel-masking studies, cluster/event
statistics ("q-statistic"), and rate plots comparing masking strategies.

## Structure

```
code/
├── config/         # config.yaml — OHDUs, pixel mask bitmasks, and input .root file lists
│                   # config_test.yaml, gif_cfg.yaml — inputs for the test/gif macros below
├── cxx/            # ROOT (C++) analysis macros
│   ├── helperFunctions.C  # shared utilities included by the macros below
│   ├── LEC.C               # "Local Event Cut" masking analysis
│   ├── multichannel.C      # multi-channel/multi-OHDU analysis (full threshold scan +
│   │                       # threshold=0.1 diagnostic pass, writes to txts/ and results/)
│   ├── binomialMask.C      # binomial pixel-mask analysis
│   ├── hduPlotter.C        # per-HDU plotting helper
│   ├── multichannel_test.C, LEC_test.C  # per-LTA/OHDU dev variants (config_test.yaml)
│   ├── run_thr0p1.C        # standalone rerun of the threshold=0.1 diagnostic pass
│   ├── scanCoverage.C      # visualizes sliding-window scan coverage for one quadrant
│   ├── q1MaskGif.C         # animates the q1 masking cut across thresholds (config: gif_cfg.yaml)
│   ├── txts/                # rate/survival tables (more_*.txt, LEC_*.txt) produced by the macros
│   ├── results/              # consolidated pdfs/, root/, and txts/ output of the macros above
│   └── gif/                  # animated gif output of q1MaskGif.C
├── plots/          # rendered outputs
│   ├── pngs/, pngs_2023/   # per-LTA/OHDU exposure maps (with nolec/wbinom/wlec variants)
│   └── qdistro/             # q-statistic distributions by exposure (0/72000/216000) and window
└── python/         # rates.py — final rate/upper-limit summary plot (rates.pdf) from the .txt tables
```

## Pipeline

1. `config/config.yaml` defines which OHDUs, pixel masks (neighbor, bleed, halo,
   crosstalk, edge, hot pixel/column, LEC, etc.), and input `.root` files to use.
2. The ROOT macros in `cxx/` (`LEC.C`, `multichannel.C`, `binomialMask.C`) process
   the hit data, apply the configured masks, and write survival/rate tables to
   `cxx/txts/` plus diagnostic PDFs and root files under `cxx/results/`.
3. `python/rates.py` reads the `more_*.txt` / `LEC_*.txt` tables and produces the
   summary rate plot (`rates.pdf`).
4. `run_thr0p1.C`, `scanCoverage.C`, and `q1MaskGif.C` are diagnostic tools built on
   top of that pipeline: rerunning just the threshold=0.1 pass, visualizing sliding-window
   scan coverage, and animating how the q1 masking cut evolves across thresholds.

## Dependencies

- [ROOT](https://root.cern) (CERN)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
- Python: `matplotlib`, `numpy`, `scipy`
