"""
compare_psf.py
==============
Multi-image PSF quality diagnostic figure.

Compares ePSF (EPSFBuilder) across 5 representative images spanning different
observing conditions, filters and seeing.  A PSFEx reference kernel is shown
for an apples-to-apples radial-profile comparison.

Layout (3 rows):
  Row 0  – 2-D PSF kernels (log scale): 5 ePSFs + PSFEx reference
  Row 1  – Radial profiles (all 6 overlaid) | 3-D surface (best image)
  Row 2  – Summary metrics table

Wing-degradation note
---------------------
With kernel_size = 6×FWHM (~27 px), the ePSF's radial profile turned upward
past r ≈ 6 px because the truncated wings were renormalised to unit sum.
The fix in build_psf.py (kernel_size = 10×FWHM, ≥51 px) matches PSFEx's
59×59 default and captures the full stellar halo.

Run from the subphot_pipe directory:
    python compare_psf.py
"""

import os
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm, PowerNorm
from mpl_toolkits.mplot3d import Axes3D           # noqa: F401
from scipy.ndimage import zoom, gaussian_filter
from astropy.io import fits
from astropy.stats import mad_std

from build_psf import (
    build_psf_from_fits,
    _fit_psf_fwhm,
)

# ─────────────────────────────────────────────────────────────────────────────
# Image set definition
# ─────────────────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))

# Each entry: (label, short_label, path, color, min_snr)
IMAGES = [
    ('ZTF26 i  240 s  2026-03-29',  'i 240 s\n03-29',
     'data/ZTF26aakjzdt/rc20260329_06_34_50_f_a_b_ZTF26aakjzdt_i_i.fits',
     '#1a7ab5', 10),   # blue  – best image

    ('ZTF26 i  240 s  2026-03-15',  'i 240 s\n03-15',
     'data/ZTF26aakjzdt/rc20260315_06_12_08_f_a_b_ZTF26aakjzdt_i_i.fits',
     '#4cb8e0', 10),   # light-blue – good, diff. seeing

    ('ZTF26 g  180 s  2026-03-29',  'g 180 s\n03-29',
     'data/ZTF26aakjzdt/rc20260329_05_13_30_f_a_b_ZTF26aakjzdt_g_g.fits',
     '#2ca02c', 10),   # green

    ('ZTF26 r  180 s  2026-03-15',  'r 180 s\n03-15',
     'data/ZTF26aakjzdt/rc20260315_06_08_19_f_a_b_ZTF26aakjzdt_r_r.fits',
     '#ff7f0e', 10),   # orange – mediocre r-band

    ('ZTF26 r  180 s  2026-03-29\n[PSFEx χ²≫1, all stars faint]',
     'r 180 s\n03-29\n(bad)',
     'data/ZTF26aakjzdt/rc20260329_05_17_21_f_a_b_ZTF26aakjzdt_r_r.fits',
     '#d62728', 0),    # red   – PSFEx-failed image (min_snr=0 fallback)
]

PSFEX_PATH  = os.path.join(BASE, 'out/proto_ref_prepsfex_33207.fits')
PSFEX_LABEL = 'PSFEx ref\nchi²=0.99  FWHM=4.18 px'
C_PSFEX     = '#9467bd'   # purple
GRAY        = '#555555'

# ─────────────────────────────────────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────────────────────────────────────

def _load_psfex_kernel(path):
    """Load PSFEx proto kernel, normalise to sum=1."""
    with fits.open(path) as hdul:
        kernel = hdul[0].data[0].astype(float)
    s = kernel.sum()
    if s > 0:
        kernel /= s
    return kernel


def _radial_profile(kernel, n_bins=50):
    """Compute azimuthally-averaged radial profile, normalised to peak."""
    ny, nx = kernel.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - cx, yy - cy)
    r_flat, k_flat = r.ravel(), kernel.ravel()
    r_max = r_flat.max()
    bins = np.linspace(0, r_max, n_bins + 1)
    rad, prof = [], []
    for i in range(n_bins):
        mask = (r_flat >= bins[i]) & (r_flat < bins[i + 1])
        if mask.sum() > 0:
            rad.append(0.5 * (bins[i] + bins[i + 1]))
            prof.append(k_flat[mask].mean())
    rad, prof = np.array(rad), np.array(prof)
    peak = prof.max()
    if peak > 0:
        prof /= peak
    return rad, prof


def _psf_panel(ax, kernel, title, color_title='black', fwhm=None, elong=None, n_stars=None):
    """Draw a 2-D PSF on a log scale with metric annotations."""
    k = kernel.copy()
    kpos = k[k > 0]
    vmin = kpos.min() if kpos.size else 1e-9
    norm = LogNorm(vmin=vmin, vmax=k.max())
    ax.imshow(k, origin='lower', cmap='inferno', norm=norm,
              interpolation='nearest', aspect='equal')
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # Kernel size label
    ny, nx = k.shape
    ax.text(0.02, 0.02, f'{nx}×{ny} px',
            transform=ax.transAxes, fontsize=6.5, color='white', va='bottom',
            bbox=dict(facecolor='black', alpha=0.35, pad=1.5, edgecolor='none'))

    # Metric strip at top
    meta = []
    if fwhm is not None:
        meta.append(f'FWHM {fwhm:.2f}')
    if elong is not None:
        meta.append(f'elong {elong:.2f}')
    if n_stars is not None:
        meta.append(f'N={n_stars}')
    if meta:
        ax.text(0.5, 1.02, '  ·  '.join(meta),
                transform=ax.transAxes, fontsize=7, ha='center', va='bottom',
                color=color_title)

    ax.set_title(title, color=color_title, fontsize=8, fontweight='bold', pad=14)
    for sp in ax.spines.values():
        sp.set_edgecolor(color_title)
        sp.set_linewidth(1.5)


def _surface3d(ax3d, kernel, color, title):
    """3-D surface plot of the PSF kernel."""
    ny, nx = kernel.shape
    yy, xx = np.mgrid[:ny, :nx]
    k = kernel / kernel.max()
    k_sm = gaussian_filter(k, sigma=0.5)
    ax3d.plot_surface(xx, yy, k_sm, cmap='inferno',
                      edgecolor='none', alpha=0.90,
                      linewidth=0, antialiased=True)
    ax3d.set_title(title, color=color, fontweight='bold', pad=2, fontsize=8)
    ax3d.set_xlabel('x [px]', fontsize=6.5, labelpad=1)
    ax3d.set_ylabel('y [px]', fontsize=6.5, labelpad=1)
    ax3d.set_zlabel('Norm flux', fontsize=6.5, labelpad=1)
    ax3d.tick_params(labelsize=6, pad=0)
    ax3d.view_init(elev=30, azim=-60)
    ax3d.set_facecolor('white')


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("Building ePSFs …")
    results = []
    for label, short_label, path, color, min_snr in IMAGES:
        abspath = os.path.join(BASE, path)
        try:
            r = build_psf_from_fits(
                abspath,
                fwhm_guess=4.5,
                threshold_sigma=4.0,
                oversampling=4,
                max_stars=40,
                min_stars=3,
                epsf_iters=3,
                epsf_smoothing_kernel='quartic',
                min_snr=min_snr,
            )
            r['label'] = label
            r['short_label'] = short_label
            r['color'] = color
            r['min_snr'] = min_snr
            r['path'] = abspath
            print(f"  {label[:40]:40s}  FWHM={r['fwhm']:.2f}  elong={r['elongation']:.2f}"
                  f"  n={r['n_stars']}  scatter={r['fwhm_scatter']:.3f}")
            results.append(r)
        except Exception as e:
            print(f"  {label[:40]:40s}  FAILED: {e}")

    psfex_kernel = _load_psfex_kernel(PSFEX_PATH)
    psfex_fx, psfex_fy, _ = _fit_psf_fwhm(psfex_kernel)
    psfex_fwhm = float(np.sqrt(psfex_fx * psfex_fy))
    print(f"  PSFEx reference                            FWHM={psfex_fwhm:.2f}  "
          f"kernel={psfex_kernel.shape[0]}×{psfex_kernel.shape[1]} px")

    n = len(results)
    if n == 0:
        raise RuntimeError("No images processed successfully.")

    print("Building figure …")

    plt.rcParams.update({
        'font.family':      'DejaVu Sans',
        'font.size':        9,
        'axes.labelsize':   9,
        'axes.titlesize':   10,
        'xtick.labelsize':  8,
        'ytick.labelsize':  8,
        'figure.facecolor': 'white',
        'axes.facecolor':   '#f8f8f8',
    })

    n_cols = n + 1   # images + PSFEx reference
    fig_width = max(18, n_cols * 3.2)

    fig = plt.figure(figsize=(fig_width, 18))
    fig.patch.set_facecolor('white')

    outer = gridspec.GridSpec(
        3, 1, figure=fig,
        top=0.93, bottom=0.04,
        left=0.05, right=0.97,
        hspace=0.45,
        height_ratios=[1.1, 1.2, 0.7],
    )

    # ─────────────────────────────────────────────────────────────────────────
    # ROW 0 – 2-D PSF kernels
    # ─────────────────────────────────────────────────────────────────────────
    gs0 = gridspec.GridSpecFromSubplotSpec(
        1, n_cols, subplot_spec=outer[0], wspace=0.3,
    )

    for ci, r in enumerate(results):
        ax = fig.add_subplot(gs0[0, ci])
        _psf_panel(ax, r['kernel'], r['short_label'], color_title=r['color'],
                   fwhm=r['fwhm'], elong=r['elongation'], n_stars=r['n_stars'])

    # PSFEx reference panel
    ax_psfex = fig.add_subplot(gs0[0, n])
    _psf_panel(ax_psfex, psfex_kernel, PSFEX_LABEL, color_title=C_PSFEX,
               fwhm=psfex_fwhm)

    # Row 0 title
    fig.text(0.5, outer[0].get_position(fig).y1 + 0.01,
             'PSF kernels (log scale)  —  51×51 px native pixels  |  PSFEx: 59×59 px',
             ha='center', va='bottom', fontsize=9.5, color=GRAY, style='italic')

    # ─────────────────────────────────────────────────────────────────────────
    # ROW 1 – Radial profiles + 3-D surface
    # ─────────────────────────────────────────────────────────────────────────
    gs1 = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[1], wspace=0.35,
        width_ratios=[1.7, 1.0],
    )

    # Radial profiles
    ax_rad = fig.add_subplot(gs1[0, 0])
    ax_rad.set_facecolor('#f8f8f8')

    for r in results:
        rad, prof = _radial_profile(r['kernel'])
        ls = '--' if r['min_snr'] == 0 else '-'
        lw = 2.2 if r['fwhm_scatter'] < 0.5 else 1.5
        label = (f"{r['short_label'].replace(chr(10),' ')}  "
                 f"FWHM={r['fwhm']:.2f} px"
                 + (f"  scatter={r['fwhm_scatter']:.2f}" if not np.isnan(r['fwhm_scatter']) else ''))
        ax_rad.semilogy(rad, np.clip(prof, 1e-5, None),
                        color=r['color'], lw=lw, ls=ls, label=label, alpha=0.9)

    # PSFEx reference
    r_p, p_p = _radial_profile(psfex_kernel)
    ax_rad.semilogy(r_p, np.clip(p_p, 1e-5, None),
                    color=C_PSFEX, lw=2.5, ls=':', zorder=5,
                    label=f'PSFEx ref  FWHM={psfex_fwhm:.2f} px  chi²=0.99')

    ax_rad.axhline(0.5, ls=':', lw=0.8, color=GRAY, alpha=0.7, label='50% level')
    ax_rad.axhline(0.1, ls=':', lw=0.6, color=GRAY, alpha=0.4)

    ax_rad.set_xlabel('Radius [pixels]', fontsize=9)
    ax_rad.set_ylabel('Normalised flux  (log scale)', fontsize=9)
    ax_rad.set_title('Radial PSF profiles  —  normalised to peak\n'
                     '(dashed = PSFEx-failed image;  thick = low scatter)',
                     fontweight='bold')
    ax_rad.set_xlim(0, 26)
    ax_rad.set_ylim(1e-4, 2)
    ax_rad.legend(fontsize=7.5, framealpha=0.9, loc='upper right')
    for sp in ['top', 'right']:
        ax_rad.spines[sp].set_visible(False)

    # Annotate the old truncation radius
    ax_rad.axvline(13, ls='--', lw=0.8, color='#aaaaaa', alpha=0.7)
    ax_rad.text(13.2, 3e-4, 'Old 6×FWHM\ntruncation\n(r=13 px)',
                fontsize=6.5, color='#aaaaaa', va='bottom')

    # 3-D surface of best image (lowest fwhm_scatter)
    best = min(results, key=lambda r: r['fwhm_scatter']
               if not np.isnan(r['fwhm_scatter']) else 99)
    ax3d = fig.add_subplot(gs1[0, 1], projection='3d')
    _surface3d(ax3d, best['kernel'], best['color'],
               f"3-D ePSF  — best image\n{best['short_label'].replace(chr(10),' ')}")

    # ─────────────────────────────────────────────────────────────────────────
    # ROW 2 – Summary metrics table
    # ─────────────────────────────────────────────────────────────────────────
    ax_tbl = fig.add_subplot(outer[2])
    ax_tbl.set_facecolor('white')
    ax_tbl.axis('off')

    # Build table
    col_labels = ['Metric'] + [r['short_label'].replace('\n', ' ') for r in results] + ['PSFEx ref']
    rows = []

    rows.append(
        ['FWHM [px]']
        + [f"{r['fwhm']:.2f}" for r in results]
        + [f'{psfex_fwhm:.2f}']
    )
    rows.append(
        ['Elongation']
        + [f"{r['elongation']:.2f}" for r in results]
        + ['—']
    )
    rows.append(
        ['N stars in ePSF']
        + [str(r['n_stars']) for r in results]
        + ['15 (est)']
    )
    rows.append(
        ['FWHM scatter [px]']
        + [(f"{r['fwhm_scatter']:.3f}" if not np.isnan(r['fwhm_scatter']) else 'NaN')
           for r in results]
        + ['—']
    )
    rows.append(
        ['Residual RMS [%]']
        + [(f"{r['residual_rms_pct']:.1f}" if not np.isnan(r['residual_rms_pct']) else 'NaN')
           for r in results]
        + ['—']
    )
    rows.append(
        ['Kernel size [px]']
        + [f"{r['kernel'].shape[0]}×{r['kernel'].shape[0]}" for r in results]
        + [f"{psfex_kernel.shape[0]}×{psfex_kernel.shape[0]}"]
    )
    rows.append(
        ['Quality']
        + [('✓ Good' if (not np.isnan(r['fwhm_scatter']) and r['fwhm_scatter'] < 0.5
                        and r['elongation'] < 1.5)
            else ('⚠ Marginal' if r['elongation'] < 2.0 else '✗ Failed'))
           for r in results]
        + ['✓ Reference']
    )

    tbl = ax_tbl.table(
        cellText=rows,
        colLabels=col_labels,
        cellLoc='center',
        loc='center',
        bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1.0, 1.6)

    n_cols_tbl = n + 2

    # Header row style
    for ci in range(n_cols_tbl):
        cell = tbl[(0, ci)]
        cell.set_facecolor('#dde3ec')
        cell.set_text_props(fontweight='bold', fontsize=7.5)

    # Metric label column
    for ri in range(1, len(rows) + 1):
        tbl[(ri, 0)].set_facecolor('#f0f0f0')
        tbl[(ri, 0)].set_text_props(fontweight='bold')

    # Image columns — colour by quality (fwhm_scatter)
    for ci, r in enumerate(results, start=1):
        sc = r['fwhm_scatter']
        if np.isnan(sc) or sc > 3.0:
            fc = '#fde8e8'  # red tint (poor)
        elif sc < 0.5:
            fc = '#e1f0fa'  # blue tint (good)
        else:
            fc = '#fff8e1'  # yellow tint (mediocre)
        for ri in range(1, len(rows) + 1):
            tbl[(ri, ci)].set_facecolor(fc)

    # PSFEx column
    for ri in range(1, len(rows) + 1):
        tbl[(ri, n + 1)].set_facecolor('#ede9f8')  # purple tint

    ax_tbl.set_title('Summary metrics  —  colour: blue=good  yellow=mediocre  red=poor/failed',
                     fontweight='bold', fontsize=9, pad=4)

    # ─────────────────────────────────────────────────────────────────────────
    # Overall title
    # ─────────────────────────────────────────────────────────────────────────
    fig.suptitle(
        'PSF Quality Assessment  —  EPSFBuilder (10×FWHM kernel) vs PSFEx\n'
        'ZTF26aakjzdt  ·  SEDM / P60  ·  Multiple epochs & filters',
        fontsize=12, fontweight='bold', y=0.96,
    )

    outpath = os.path.join(BASE, 'psf_comparison.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"✓ Saved: {outpath}")
    return outpath


if __name__ == '__main__':
    main()
