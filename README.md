# sensei2026

Analysis code for SENSEI Skipper-CCD data: pixel-masking studies, cluster/event
statistics ("q-statistic"), and rate plots comparing masking strategies.

## Structure

```
code/
├── config/         # config.yaml — OHDUs, pixel mask bitmasks, and input .root file lists
├── cxx/            # ROOT (C++) analysis macros
│   ├── helperFunctions.C  # shared utilities included by the macros below
│   ├── LEC.C               # "Local Event Cut" masking analysis
│   ├── multichannel.C      # multi-channel/multi-OHDU analysis
│   ├── binomialMask.C      # binomial pixel-mask analysis
│   ├── txts/                # rate/survival tables produced by the macros, by run/window (r8-r30)
│   └── pdfs/, EXPO216000/, recent/  # q-statistic and mask plots (PDF) per window/exposure
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
   `cxx/txts/` plus diagnostic PDFs.
3. `python/rates.py` reads the `more_*.txt` / `LEC_*.txt` tables and produces the
   summary rate plot (`rates.pdf`).

## Dependencies

- [ROOT](https://root.cern) (CERN)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
- Python: `matplotlib`, `numpy`, `scipy`
