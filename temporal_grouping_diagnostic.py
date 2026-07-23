"""
Diagnostic tools for understanding why temporal grouping is failing.
"""

import pandas as pd
import numpy as np


def diagnose_temporal_overlap(df, instruments=['NOT/ALFOSC', 'LT/IOO'], time_window=1.0):
    """
    Show why images aren't being grouped together.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with obs_mjd, instrument columns
    instruments : list
        Instruments to check
    time_window : float
        Time window in days
    """
    print(f"Diagnostic: Why aren't {instruments[0]} and {instruments[1]} being grouped?\n")

    inst_dfs = {}
    for inst in instruments:
        inst_dfs[inst] = df[df['instrument'] == inst].sort_values('obs_mjd')
        if len(inst_dfs[inst]) > 0:
            print(f"{inst}:")
            print(f"  Count: {len(inst_dfs[inst])}")
            print(f"  MJD range: {inst_dfs[inst]['obs_mjd'].min():.3f} - {inst_dfs[inst]['obs_mjd'].max():.3f}")
            print(f"  Span: {inst_dfs[inst]['obs_mjd'].max() - inst_dfs[inst]['obs_mjd'].min():.3f} days")
            print()

    # Check if there's any temporal overlap
    if len(instruments) == 2:
        inst1, inst2 = instruments
        if len(inst_dfs[inst1]) > 0 and len(inst_dfs[inst2]) > 0:
            min1, max1 = inst_dfs[inst1]['obs_mjd'].min(), inst_dfs[inst1]['obs_mjd'].max()
            min2, max2 = inst_dfs[inst2]['obs_mjd'].min(), inst_dfs[inst2]['obs_mjd'].max()

            print("Temporal overlap check:")
            if max1 < min2 - time_window:
                print(f"  ✗ {inst1} observations end before {inst2} starts")
                print(f"    {inst1} max: {max1:.3f}, {inst2} min: {min2:.3f}, gap: {min2 - max1:.3f} days")
            elif max2 < min1 - time_window:
                print(f"  ✗ {inst2} observations end before {inst1} starts")
                print(f"    {inst2} max: {max2:.3f}, {inst1} min: {min1:.3f}, gap: {min1 - max2:.3f} days")
            elif max1 >= min2 - time_window and max2 >= min1 - time_window:
                print(f"  ✓ Observations overlap within {time_window} day window")

                # Show which observations are close
                print(f"\n  Closest observations:")
                for mjd1, obj1 in zip(inst_dfs[inst1]['obs_mjd'], inst_dfs[inst1]['object'])[:3]:
                    diffs = np.abs(inst_dfs[inst2]['obs_mjd'] - mjd1)
                    closest_idx = diffs.argmin()
                    mjd2 = inst_dfs[inst2].iloc[closest_idx]['obs_mjd']
                    obj2 = inst_dfs[inst2].iloc[closest_idx]['object']
                    diff = mjd2 - mjd1
                    print(f"    {inst1} @ {mjd1:.3f} ({obj1}) <-> {inst2} @ {mjd2:.3f} ({obj2}): {abs(diff):.4f} days")


def show_mjd_gaps(df, instruments=['NOT/ALFOSC', 'LT/IOO']):
    """
    Show MJD distribution for each instrument.
    """
    print("\nMJD gaps and distribution:\n")

    for inst in instruments:
        inst_df = df[df['instrument'] == inst].sort_values('obs_mjd')
        if len(inst_df) > 1:
            mjds = inst_df['obs_mjd'].values
            gaps = np.diff(mjds)

            print(f"{inst}:")
            print(f"  Observations: {len(mjds)}")
            if len(gaps) > 0:
                print(f"  Gap between obs: min={gaps.min():.4f}, max={gaps.max():.4f}, mean={gaps.mean():.4f} days")
            print()


def suggest_time_window(df, instruments=['NOT/ALFOSC', 'LT/IOO']):
    """
    Suggest an appropriate time window for grouping.
    """
    print("\nSuggested time windows:\n")

    for inst in instruments:
        inst_df = df[df['instrument'] == inst].sort_values('obs_mjd')
        if len(inst_df) > 1:
            mjds = inst_df['obs_mjd'].values
            gaps = np.diff(mjds)

            # Time when observations are taken (e.g., within a night)
            within_night_gap = gaps[gaps < 0.5].max() if np.any(gaps < 0.5) else None
            # Time between observing runs
            between_run_gap = gaps[gaps > 0.5].min() if np.any(gaps > 0.5) else None

            print(f"{inst}:")
            if within_night_gap:
                print(f"  Within-night observations separated by up to {within_night_gap:.4f} days")
            if between_run_gap:
                print(f"  Observing runs separated by at least {between_run_gap:.4f} days")
            print()
