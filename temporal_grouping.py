"""
Group multi-telescope observations by temporal proximity.

Finds images from different instruments that are taken within a specified
time window (e.g., same night, within 1 day, etc.) for multi-filter color composites.

Usage:
    from temporal_grouping import find_temporal_groups, extract_object_name

    df['object'] = df['filename'].apply(extract_object_name)
    groups = find_temporal_groups(df, time_window=1.0, instruments=['NOT/ALFOSC', 'LT/IOO'])
"""

import re
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def extract_object_name(filepath: str) -> str:
    """
    Extract object name from file path.

    Handles patterns like:
    - /path/to/ZTF26aakjzdt_LT/image.fits → ZTF26aakjzdt
    - /path/to/ZTF25aceeneu/image.fits → ZTF25aceeneu
    - /path/to/NOT_20260423-24/image.fits → NOT_20260423-24

    Parameters
    ----------
    filepath : str
        Path to FITS file

    Returns
    -------
    object_name : str
        Extracted object identifier
    """
    path = Path(filepath)

    # Try to find ZTF name (most reliable)
    match = re.search(r'(ZTF\d+[a-z]+)', filepath)
    if match:
        return match.group(1)

    # Try directory name above file
    parent = path.parent.name
    if parent and parent != '.':
        # Remove instrument suffixes
        obj = re.sub(r'_(LT|NOT|TJO)$', '', parent)
        # Remove dates
        if 'NOT_' not in obj:
            return obj

    # Fall back to parent directory
    return path.parent.name


def find_temporal_groups(df: pd.DataFrame,
                        time_window: float = 1.0,
                        instruments: Optional[List[str]] = None,
                        require_filters: Optional[List[str]] = None) -> Dict[str, List[pd.DataFrame]]:
    """
    Group observations by object and temporal proximity.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns: 'filename', 'obs_mjd', 'instrument', 'filter'
        Optionally: 'object' (if not present, will be extracted)
    time_window : float
        Time window in days to consider images as simultaneous
    instruments : list, optional
        Only consider these instruments (e.g., ['NOT/ALFOSC', 'LT/IOO'])
    require_filters : list, optional
        Only consider these filters (e.g., ['g', 'r', 'i'])

    Returns
    -------
    groups : dict
        {object_name: [list of DataFrames for valid groups]}
        Each DataFrame is a group of images within time_window with
        images from different instruments
    """
    df = df.copy()

    # Add object column if not present
    if 'object' not in df.columns:
        df['object'] = df['filename'].apply(extract_object_name)

    # Filter by instrument if specified
    if instruments:
        df = df[df['instrument'].isin(instruments)]

    # Filter by filters if specified
    if require_filters:
        df = df[df['filter'].isin(require_filters)]

    if len(df) == 0:
        print("No data matching filters")
        return {}

    groups = {}

    # Group by object
    for obj_name, obj_df in df.groupby('object'):
        obj_groups = []

        # Sort by MJD
        obj_df = obj_df.sort_values('obs_mjd').reset_index(drop=True)

        # Find temporal groups
        used = set()
        for i, (idx, row) in enumerate(obj_df.iterrows()):
            if idx in used:
                continue

            # Find all images within time_window of this one
            mjd_center = row['obs_mjd']
            mask = np.abs(obj_df['obs_mjd'] - mjd_center) <= time_window
            group = obj_df[mask].copy()

            # Only keep groups with multiple instruments
            if len(group['instrument'].unique()) > 1:
                obj_groups.append(group)
                used.update(group.index)

        if obj_groups:
            groups[obj_name] = obj_groups

    return groups


def summarize_groups(groups: Dict[str, List[pd.DataFrame]]) -> pd.DataFrame:
    """
    Summarize temporal groups as a table.

    Parameters
    ----------
    groups : dict
        Output from find_temporal_groups()

    Returns
    -------
    summary : pd.DataFrame
        Table with columns: object, group_id, num_images, filters, instruments,
        mjd_min, mjd_max, mjd_span
    """
    rows = []

    for obj_name, obj_groups in groups.items():
        for group_id, group_df in enumerate(obj_groups):
            rows.append({
                'object': obj_name,
                'group_id': group_id,
                'num_images': len(group_df),
                'filters': ','.join(sorted(group_df['filter'].unique())),
                'instruments': ','.join(sorted(group_df['instrument'].unique())),
                'mjd_min': group_df['obs_mjd'].min(),
                'mjd_max': group_df['obs_mjd'].max(),
                'mjd_span_days': group_df['obs_mjd'].max() - group_df['obs_mjd'].min(),
            })

    return pd.DataFrame(rows)


def print_groups(groups: Dict[str, List[pd.DataFrame]], verbose: bool = True) -> None:
    """
    Print a human-readable summary of temporal groups.

    Parameters
    ----------
    groups : dict
        Output from find_temporal_groups()
    verbose : bool
        Print detailed info about each image in each group
    """
    if not groups:
        print("No temporal groups found")
        return

    print(f"Found {sum(len(v) for v in groups.values())} temporal groups\n")

    for obj_name in sorted(groups.keys()):
        obj_groups = groups[obj_name]
        print(f"{'='*70}")
        print(f"Object: {obj_name}  ({len(obj_groups)} group{'s' if len(obj_groups) > 1 else ''})")
        print(f"{'='*70}")

        for group_id, group_df in enumerate(obj_groups):
            print(f"\n  Group {group_id + 1}:")
            print(f"    Time span: {group_df['obs_mjd'].max() - group_df['obs_mjd'].min():.4f} days")
            print(f"    Filters: {','.join(sorted(group_df['filter'].unique()))}")
            print(f"    Instruments: {','.join(sorted(group_df['instrument'].unique()))}")
            print(f"    Images: {len(group_df)}")

            if verbose:
                print(f"\n    Details:")
                for _, img in group_df.iterrows():
                    inst_short = img['instrument'].split('/')[-1]
                    print(f"      {img['obs_mjd']:.4f}  {img['filter']}  {inst_short:6s}  "
                          f"{Path(img['filename']).name}")


def select_group(groups: Dict[str, List[pd.DataFrame]],
                object_name: str,
                group_id: int = 0) -> pd.DataFrame:
    """
    Select a single temporal group.

    Parameters
    ----------
    groups : dict
        Output from find_temporal_groups()
    object_name : str
        Which object
    group_id : int
        Which group (0-indexed)

    Returns
    -------
    group_df : pd.DataFrame
        Selected group as DataFrame
    """
    if object_name not in groups:
        raise ValueError(f"Object '{object_name}' not in groups. "
                        f"Available: {list(groups.keys())}")

    if group_id >= len(groups[object_name]):
        raise ValueError(f"Group {group_id} not found for {object_name}. "
                        f"Only {len(groups[object_name])} groups available")

    return groups[object_name][group_id]


def plot_temporal_distribution(df: pd.DataFrame,
                              figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    Plot temporal distribution of observations.

    Shows which images are taken close in time across instruments.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'obs_mjd', 'instrument', 'filter', 'object'
    figsize : tuple
        Figure size

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    if 'object' not in df.columns:
        df = df.copy()
        df['object'] = df['filename'].apply(extract_object_name)

    fig, axes = plt.subplots(len(df['object'].unique()), 1,
                            figsize=figsize, sharex=True)
    if len(df['object'].unique()) == 1:
        axes = [axes]

    for ax, obj_name in enumerate(sorted(df['object'].unique())):
        obj_df = df[df['object'] == obj_name]

        # Plot each instrument with different color
        instruments = sorted(obj_df['instrument'].unique())
        colors = plt.cm.Set1(np.linspace(0, 1, len(instruments)))
        markers = {'g': 'o', 'r': 's', 'i': '^', 'z': 'D', 'u': 'v', 'y': 'p'}

        for inst, color in zip(instruments, colors):
            inst_df = obj_df[obj_df['instrument'] == inst]
            for filt in inst_df['filter'].unique():
                filt_df = inst_df[inst_df['filter'] == filt]
                marker = markers.get(filt, 'o')

                axes[ax].scatter(filt_df['obs_mjd'], [inst] * len(filt_df),
                               color=color, marker=marker, s=100, alpha=0.7,
                               label=f'{filt}')

        axes[ax].set_ylabel(obj_name, fontsize=10)
        axes[ax].set_yticks(instruments)
        axes[ax].grid(True, alpha=0.3)

    axes[-1].set_xlabel('MJD', fontsize=11)
    fig.suptitle('Temporal Distribution of Observations', fontsize=13, fontweight='bold')
    plt.tight_layout()

    return fig


# Quick utility to show what filters are available per group
def filter_coverage(group_df: pd.DataFrame) -> Dict[str, List[str]]:
    """
    Show which filters are available for each instrument in a group.

    Parameters
    ----------
    group_df : pd.DataFrame
        Single temporal group

    Returns
    -------
    coverage : dict
        {instrument: [list of filters]}
    """
    coverage = {}
    for inst in group_df['instrument'].unique():
        filters = sorted(group_df[group_df['instrument'] == inst]['filter'].unique())
        coverage[inst] = filters
    return coverage
