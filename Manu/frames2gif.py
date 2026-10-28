#!/usr/bin/env python3
"""Assemble MaskAnimation.C's frame_NN.png files into an animated GIF.

usage: frames2gif.py <frames_dir> [out.gif] [--width 1200]
The first three frames (raw image, mask, grid) and the last one (result
maps) are held longer than the window-by-window frames.
"""
import glob, os, sys
from PIL import Image

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    width = 1200
    if '--width' in sys.argv:
        width = int(sys.argv[sys.argv.index('--width') + 1])
        args = [a for a in args if a != str(width)]
    if not args:
        sys.exit(__doc__)
    d = args[0]
    out = args[1] if len(args) > 1 else os.path.join(d, 'mask_walkthrough.gif')

    files = sorted(glob.glob(os.path.join(d, 'frame_*.png')))
    if not files:
        sys.exit('no frame_*.png in ' + d)

    frames = []
    for f in files:
        im = Image.open(f).convert('RGB')
        h = round(im.height * width / im.width)
        im = im.resize((width, h), Image.LANCZOS)
        frames.append(im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))

    n = len(frames)
    durations = [3500 if (i < 3 or i == n - 1) else 2200 for i in range(n)]
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=durations,
                   loop=0, optimize=True, disposal=1)
    print('wrote', out, f'({n} frames, {os.path.getsize(out)/1e6:.1f} MB)')

if __name__ == '__main__':
    main()
