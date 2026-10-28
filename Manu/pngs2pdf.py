#!/usr/bin/env python3
"""Combine PNG images into one PDF without quality loss.

Each image becomes one page, stored with lossless Flate compression (Pillow's
own PDF writer stores colour images as JPEG, which smears single pixels and
blurs text).

usage: pngs2pdf.py out.pdf page1.png page2.png ... [--dpi 300]
"""
import sys
import zlib
from PIL import Image


def write_pdf(out, files, dpi):
    objs = []                                   # body of object i+1

    def add(body):
        objs.append(body)
        return len(objs)

    catalog = add(None)
    pages = add(None)
    kids = []

    for f in files:
        im = Image.open(f).convert('RGB')
        w, h = im.size
        data = zlib.compress(im.tobytes(), 9)
        img = add(b'<< /Type /XObject /Subtype /Image /Width %d /Height %d /ColorSpace /DeviceRGB '
                  b'/BitsPerComponent 8 /Filter /FlateDecode /Length %d >>\nstream\n' % (w, h, len(data))
                  + data + b'\nendstream')
        pw, ph = w * 72.0 / dpi, h * 72.0 / dpi
        content = b'q %.4f 0 0 %.4f 0 0 cm /Im0 Do Q' % (pw, ph)
        cont = add(b'<< /Length %d >>\nstream\n' % len(content) + content + b'\nendstream')
        page = add(b'<< /Type /Page /Parent %d 0 R /MediaBox [0 0 %.4f %.4f] '
                   b'/Resources << /XObject << /Im0 %d 0 R >> >> /Contents %d 0 R >>' % (pages, pw, ph, img, cont))
        kids.append(page)

    objs[catalog - 1] = b'<< /Type /Catalog /Pages %d 0 R >>' % pages
    objs[pages - 1] = (b'<< /Type /Pages /Count %d /Kids [' % len(kids)
                       + b' '.join(b'%d 0 R' % k for k in kids) + b'] >>')

    with open(out, 'wb') as fh:
        fh.write(b'%PDF-1.4\n%\xe2\xe3\xcf\xd3\n')
        offsets = []
        for i, body in enumerate(objs, 1):
            offsets.append(fh.tell())
            fh.write(b'%d 0 obj\n' % i + body + b'\nendobj\n')
        xref = fh.tell()
        fh.write(b'xref\n0 %d\n0000000000 65535 f \n' % (len(objs) + 1))
        for off in offsets:
            fh.write(b'%010d 00000 n \n' % off)
        fh.write(b'trailer\n<< /Size %d /Root %d 0 R >>\nstartxref\n%d\n%%%%EOF\n' % (len(objs) + 1, catalog, xref))


def main():
    args = sys.argv[1:]
    dpi = 300
    if '--dpi' in args:
        i = args.index('--dpi')
        dpi = float(args[i + 1])
        del args[i:i + 2]
    if len(args) < 2:
        sys.exit(__doc__)
    write_pdf(args[0], args[1:], dpi)
    print('wrote', args[0], f'({len(args) - 1} pages)')


if __name__ == '__main__':
    main()
