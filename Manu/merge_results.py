#!/usr/bin/env python3
"""Merge the r = 30 result PDFs into one document per dataset and quad
selection, written to summary/<quads>/<dataset name>.pdf
(e.g. summary/goodquads/Old.pdf).

Each page gets a title strip on top (study name, then dataset, quad
selection and exposure), and each document gets a bookmark tree
exposure > study. Pages stay vector (Ghostscript).

usage: merge_results.py [output_dir]     (run from the project folder)
"""
import os
import re
import subprocess
import sys
import tempfile

RADIUS = 30
BANNER = 36  # points added on top of each page for the title strip

# (folder, label, output file name)
DATASETS = [
    ('.',                                'Original samples: hits_corr_proc, 2026-03-17 run 8', 'Old'),
    ('results/hits_skp_2025-09-28',      'hits_skp, 2025-09-28 (runs 28, 30, 32, 34)',         'Sept'),
    ('results/hits_skp_2026-03-17_run8', 'hits_skp, 2026-03-17 (run 8)',                       'March'),
    ('results/hits_skp_2026-07-03',      'hits_skp, 2026-07-03 (runs 5, 7)',                   'July'),
    ('results/hits_skp_all',             'hits_skp, all campaigns combined',                   'Inclusive'),
]

CATEGORIES = [('goodquads', 'good quads'), ('notbadquads', 'not-bad quads')]

STUDIES = [
    ('nocut_distributions',   'All windows (no cut): N, k and q distributions'),
    ('k2eq1_selection',       'Windows with k2 = 1: N, k and q distributions'),
    ('k2eq2_selection',       'Windows with k2 = 2: N, k and q distributions'),
    ('k3pe1_selection',       'Windows with k3+ = 1: N, k and q distributions'),
    ('qdist_nocut',           'q = k/N distributions, no cut'),
    ('qdist_k2lt2',           'q = k/N distributions, k2 < 2'),
    ('qdist_compare',         'q = k/N distributions: no cut vs k2 < 2'),
    ('qcut_scan_nocut',       'q-cut scan vs LEC mask, no cut'),
    ('qcut_scan_k2lt2',       'q-cut scan vs LEC mask, k2 < 2'),
    ('k1_independence_k2ge1', 'Independence test: k1 of windows with k2 >= 1'),
    ('k1_threecases',         'k1 distributions and P(k2 >= 1) vs k1'),
    ('k1_independence_k2fake','Closure test (images with no 2e hit, fake 2e hits): k1 independence'),
    ('k1_threecases_k2fake',  'Closure test: k1 distributions and P(fake k2 >= 1) vs k1'),
    ('k1_closure_calibration','Closure test: p-values and pulls over 20 toys'),
    ('k1_independence_k3pe1', 'Independence test: k1 of windows with k3+ = 1'),
    ('k1_threecases_k3pe1',   'k1 distributions and P(k3+ = 1) vs k1'),
    ('k2_independence_k3pe1', 'Independence test: k2 of windows with k3+ = 1'),
    ('k2_threecases_k3pe1',   'k2 distributions and P(k3+ = 1) vs k2'),
]


def ps_string(s):
    return '(' + s.replace('\\', '\\\\').replace('(', '\\(').replace(')', '\\)') + ')'


def page_size(path):
    out = subprocess.run(['pdfinfo', path], capture_output=True, text=True, check=True).stdout
    m = re.search(r'Page size:\s+([\d.]+) x ([\d.]+) pts', out)
    return float(m.group(1)), float(m.group(2))


def stamp(src, dst, title, context):
    """Copy src onto a page BANNER points taller, with a title strip on top."""
    w, h = page_size(src)
    H = h + BANNER
    ps = f'''
<< /EndPage {{
    exch pop 2 lt {{
      gsave initclip
      0.953 0.949 0.933 setrgbcolor 0 {h} {w} {BANNER} rectfill
      0.765 0.761 0.718 setrgbcolor 0.5 setlinewidth 0 {h} moveto {w} {h} lineto stroke
      0.043 0.043 0.043 setrgbcolor
      /Helvetica-Bold findfont 11 scalefont setfont
      10 {h + BANNER - 15} moveto {ps_string(title)} show
      0.322 0.318 0.306 setrgbcolor
      /Helvetica findfont 8.5 scalefont setfont
      10 {h + 8} moveto {ps_string(context)} show
      grestore true
    }} {{ false }} ifelse
  }} >> setpagedevice
'''
    subprocess.run(['gs', '-q', '-o', dst, '-sDEVICE=pdfwrite',
                    f'-dDEVICEWIDTHPOINTS={w:.2f}', f'-dDEVICEHEIGHTPOINTS={H:.2f}',
                    '-dFIXEDMEDIA', '-c', ps, '-f', src], check=True)


def merge_one(out, dlabel, clabel, groups):
    """groups: [(exposure label, [(study title, pdf path)])] -> one PDF."""
    with tempfile.TemporaryDirectory() as tmp:
        stamped = []
        marks = [f'[ /Title {ps_string(f"{dlabel}, {clabel}, r = {RADIUS}")} /DOCINFO pdfmark']
        page = 1
        for glabel, pages in groups:
            marks.append(f'[ /Count {len(pages)} /Page {page} /Title {ps_string(glabel)} /OUT pdfmark')
            for stitle, f in pages:
                dst = os.path.join(tmp, f'p{page:04d}.pdf')
                stamp(f, dst, stitle, f'{dlabel}   |   {clabel}, {glabel}   |   r = {RADIUS}')
                stamped.append(dst)
                marks.append(f'[ /Page {page} /Title {ps_string(stitle)} /OUT pdfmark')
                page += 1

        mfile = os.path.join(tmp, 'marks.ps')
        with open(mfile, 'w') as fh:
            fh.write('\n'.join(marks) + '\n')

        subprocess.run(['gs', '-q', '-o', out, '-sDEVICE=pdfwrite'] + stamped + [mfile], check=True)

    return page - 1


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else 'summary'
    os.makedirs(outdir, exist_ok=True)

    missing = []
    for ddir, dlabel, dstem in DATASETS:
        for cat, clabel in CATEGORIES:
            cdir = os.path.join(ddir, cat)
            if not os.path.isdir(cdir):
                continue
            exps = sorted(int(e[1:]) for e in os.listdir(cdir) if re.fullmatch(r'e\d+', e))
            groups = []
            for exp in exps:
                rdir = os.path.join(cdir, f'e{exp}', f'r{RADIUS}')
                tag = f'{cat}_e{exp}_r{RADIUS}'
                pages = []
                for stem, stitle in STUDIES:
                    f = os.path.join(rdir, f'{stem}_{tag}.pdf')
                    if os.path.exists(f):
                        pages.append((stitle, f))
                    else:
                        missing.append(f)
                if pages:
                    groups.append((f'exposure {exp}', pages))
            if not groups:
                continue
            os.makedirs(os.path.join(outdir, cat), exist_ok=True)
            out = os.path.join(outdir, cat, f'{dstem}.pdf')
            n = merge_one(out, dlabel, clabel, groups)
            print(f'wrote {out}: {n} pages')

    for f in missing:
        print('missing, skipped:', f)


if __name__ == '__main__':
    main()
