import argparse
import os
import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle

plt.style.use('dark_background')

NX, NY = 3200, 520

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(SCRIPT_DIR, '..', 'configs', 'qapplet_cfg.yaml')


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load(data_dir):
    image    = np.loadtxt(f"{data_dir}/image.txt")
    premask  = np.loadtxt(f"{data_dir}/premask.txt")
    windows  = pd.read_csv(f"{data_dir}/windows.csv")
    return image, premask, windows


def plot_events(ax, event_ys, event_xs, event_vals, vmax, zorder=0.5):
    """Mark individual nonzero-charge pixels (the 'events' the q1 statistic
    is built from) as small markers, sized/colored by electron count, since
    a single pixel is otherwise nearly invisible against the full image."""
    if len(event_vals) == 0:
        return
    ax.scatter(event_xs + 0.5, event_ys + 0.5,
               s=15 + 8 * (event_vals - 1), c=event_vals,
               cmap='cool', vmin=1, vmax=vmax,
               edgecolors='black', linewidths=0.3, zorder=zorder)


def build(data_dir, threshold=None, bins=60, qmax=None):
    image, premask, windows = load(data_dir)
    radius = windows['radius'].iloc[0]
    event_ys, event_xs = np.nonzero(image)
    event_vals = image[event_ys, event_xs]
    event_vmax = max(event_vals.max(), 1) if len(event_vals) else 1

    if qmax is None:
        qmax = np.quantile(windows['q1'], 0.995)
        if qmax <= 0:
            qmax = windows['q1'].max() if windows['q1'].max() > 0 else 1e-3

    fig = plt.figure(figsize=(19, 9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1.3], height_ratios=[2, 1],
                            hspace=0.35, wspace=0.3)
    ax_img  = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[1, 0])
    ax_zoom = fig.add_subplot(gs[:, 1])

    ax_img.imshow(image, origin='lower', extent=(0, NX, 0, NY),
                  cmap='inferno', vmin=0, vmax=5, aspect='auto', zorder=0)

    pre_rgba = np.zeros((*premask.shape, 4))
    pre_rgba[..., 0] = 1.0
    pre_rgba[..., 3] = premask * 0.6
    ax_img.imshow(pre_rgba, origin='lower', extent=(0, NX, 0, NY),
                  aspect='auto', zorder=1)

    plot_events(ax_img, event_ys, event_xs, event_vals, event_vmax)

    for _, w in windows[windows['masked'] == 1].iterrows():
        ax_img.add_patch(Circle((w.xcenter, w.ycenter), radius,
                                 facecolor='violet', edgecolor='none',
                                 alpha=0.5, zorder=2))

    ax_img.set_xlim(0, NX)
    ax_img.set_ylim(0, NY)
    ax_img.set_xlabel('x')
    ax_img.set_ylabel('y')
    img_default_title = 'Masked windows (q1 above config threshold)'
    ax_img.set_title(img_default_title)

    counts, bin_edges, _ = ax_hist.hist(
        windows['q1'], bins=bins, range=(0, qmax),
        color='white', histtype='step', log=True)
    ax_hist.set_xlabel(r'$q_1 = k_1/N$')
    ax_hist.set_ylabel('windows')
    hist_default_title = 'Click the q1 axis to paint windows to the right of that value'
    ax_hist.set_title(hist_default_title)

    if threshold is not None and threshold <= qmax:
        ax_hist.axvline(threshold, color='violet', lw=2,
                         label=f'mask threshold = {threshold}')
        ax_hist.legend(loc='upper right')

    # Fix the view to (0, qmax) and stop later axvlines (in on_click) from
    # re-expanding it via autoscale.
    ax_hist.set_xlim(0, qmax)
    ax_hist.set_autoscalex_on(False)

    zoom_default_title = 'Click a window in the image to zoom in'
    ax_zoom.set_title(zoom_default_title)
    ax_zoom.set_xticks([])
    ax_zoom.set_yticks([])
    ax_zoom.text(0.5, 0.5, 'no window selected', color='gray',
                 ha='center', va='center', transform=ax_zoom.transAxes)

    state = {'img_patches': [], 'hist_span': None, 'zoom_marker': None}

    def clear_highlights():
        for p in state['img_patches']:
            p.remove()
        state['img_patches'] = []
        if state['hist_span'] is not None:
            state['hist_span'].remove()
            state['hist_span'] = None
        ax_img.set_title(img_default_title)
        ax_hist.set_title(hist_default_title)

    def paint_windows(rows):
        for _, w in rows.iterrows():
            patch = Circle((w.xcenter, w.ycenter), radius,
                            facecolor='violet', edgecolor='none',
                            alpha=0.5, zorder=3)
            ax_img.add_patch(patch)
            state['img_patches'].append(patch)

    def nearest_window(x, y):
        dist2 = (windows['xcenter'] - x) ** 2 + (windows['ycenter'] - y) ** 2
        inside = windows[dist2 <= radius ** 2]
        rows = inside if len(inside) else windows.loc[[dist2.idxmin()]]
        d2 = (rows['xcenter'] - x) ** 2 + (rows['ycenter'] - y) ** 2
        return rows.loc[d2.idxmin()]

    def show_zoom(row):
        ax_zoom.cla()
        ax_zoom.set_xticks([])
        ax_zoom.set_yticks([])
        pad = radius
        x0, x1 = max(0, row.xcenter - 2 * radius - pad), min(NX, row.xcenter + 2 * radius + pad)
        y0, y1 = max(0, row.ycenter - 2 * radius - pad), min(NY, row.ycenter + 2 * radius + pad)

        ax_zoom.imshow(image, origin='lower', extent=(0, NX, 0, NY),
                       cmap='inferno', vmin=0, vmax=5, aspect='equal', zorder=0)
        ax_zoom.imshow(pre_rgba, origin='lower', extent=(0, NX, 0, NY),
                       aspect='equal', zorder=1)

        in_view = ((event_xs + 1 >= x0) & (event_xs <= x1) &
                   (event_ys + 1 >= y0) & (event_ys <= y1))
        plot_events(ax_zoom, event_ys[in_view], event_xs[in_view], event_vals[in_view], event_vmax)

        nearby = windows[windows['masked'] == 1]
        nearby = nearby[(nearby['xcenter'] + radius >= x0) & (nearby['xcenter'] - radius <= x1) &
                         (nearby['ycenter'] + radius >= y0) & (nearby['ycenter'] - radius <= y1)]
        for _, w in nearby.iterrows():
            ax_zoom.add_patch(Circle((w.xcenter, w.ycenter), radius,
                                      facecolor='violet', edgecolor='none',
                                      alpha=0.5, zorder=2))

        ax_zoom.add_patch(Circle((row.xcenter, row.ycenter), radius,
                                  facecolor='none', edgecolor='white',
                                  linewidth=2, zorder=3))
        ax_zoom.set_xlim(x0, x1)
        ax_zoom.set_ylim(y0, y1)
        ax_zoom.set_title(
            f"window @ ({row.xcenter:.0f},{row.ycenter:.0f})\n"
            f"q1={row.q1:.4f}  k1={int(row.k1)}  k2={int(row.k2)}  "
            f"k3+={int(row.k3plus)}  N={int(row.trials)}  "
            f"masked={'yes' if row.masked else 'no'}")

    def on_click(event):
        if event.xdata is None:
            return

        if event.inaxes is ax_hist:
            clear_highlights()
            q = event.xdata
            state['hist_span'] = ax_hist.axvspan(q, qmax, color='violet', alpha=0.3)
            rows = windows[windows['q1'] >= q]
            paint_windows(rows)
            ax_hist.set_title(f"q1 >= {q:.4f}: {len(rows)} windows painted")
            fig.canvas.draw_idle()

        elif event.inaxes is ax_img:
            row = nearest_window(event.xdata, event.ydata)
            if state['zoom_marker'] is not None:
                state['zoom_marker'].remove()
            marker = Circle((row.xcenter, row.ycenter), radius,
                             facecolor='none', edgecolor='white',
                             linewidth=2, zorder=4)
            ax_img.add_patch(marker)
            state['zoom_marker'] = marker
            show_zoom(row)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('button_press_event', on_click)
    return fig, (ax_img, ax_hist, ax_zoom)


NON_INTERACTIVE_BACKENDS = {'agg', 'pdf', 'ps', 'svg', 'cairo', 'template', 'pgf'}


def ensure_interactive_backend():
    """Clicks are silently no-ops if matplotlib picked a non-interactive
    backend (commonly 'agg' when no GUI toolkit / display is detected).
    Try to switch to one that actually processes mouse events."""
    current = plt.get_backend()
    if current.lower() not in NON_INTERACTIVE_BACKENDS:
        return current
    for candidate in ('QtAgg', 'TkAgg', 'GTK4Agg', 'GTK3Agg', 'MacOSX'):
        try:
            plt.switch_backend(candidate)
            return candidate
        except Exception:
            continue
    return None


if __name__ == '__main__':
    cfg = load_config()
    default_data_dir = os.path.normpath(os.path.join(SCRIPT_DIR, '..', 'cxx', cfg['output_dir']))

    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', default=default_data_dir)
    parser.add_argument('--threshold', type=float, default=cfg['threshold'])
    parser.add_argument('--bins', type=int, default=60)
    parser.add_argument('--qmax', type=float, default=None,
                         help='Upper edge of the q1 axis (default: 99.5th percentile of the data)')
    args = parser.parse_args()

    backend = ensure_interactive_backend()
    if backend is None:
        raise SystemExit(
            "No interactive matplotlib backend is available (tried Qt/Tk/GTK/MacOSX), "
            f"so it fell back to '{plt.get_backend()}', which cannot receive mouse clicks.\n"
            "Install a GUI toolkit binding, e.g.:\n"
            "  sudo apt install python3-tk        # for TkAgg\n"
            "  pip install PyQt5                  # for QtAgg\n"
            "then re-run this script.")
    print(f"[qApplet] using matplotlib backend: {backend}")

    build(args.data_dir, threshold=args.threshold, bins=args.bins, qmax=args.qmax)
    plt.show()
