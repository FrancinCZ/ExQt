"""
partitioning_plots.py — Dedicated module for Partitioning Coefficient (K_part) and Condensate Classification Analysis.
"""

from __future__ import annotations

import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


CLASS_COLOR_PALETTE = {
    "Unconstrained Droplet (Spherical LLPS)": "#2ca02c",
    "Wetted / Chromatin-Bound Condensate": "#1f77b4",
    "Hollow / Shell-Dominated": "#ff7f0e",
    "Complex / Multi-Phase Aggregate": "#d62728",
    "Punctate / Small Cluster": "#7f7f7f",
    "Unclassified": "#9467bd",
}


def plot_size_vs_partitioning(
    df: pd.DataFrame,
    output_path: Path | str | None = None,
    ax: plt.Axes | None = None,
    sample_name: str | None = None,
) -> plt.Axes:
    """Plot Condensate Size (Volume in um3) vs. Partition Coefficient (K_part)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)

    #Filter finite positive values
    size_col = "volume_bio_um3" if "volume_bio_um3" in df.columns else "shape_metric_bio"
    if size_col not in df.columns or "partition_coefficient" not in df.columns:
        if standalone:
            plt.close(fig)
        return ax

    valid = df[
        np.isfinite(df[size_col])
        & np.isfinite(df["partition_coefficient"])
        & (df[size_col] > 0)
        & (df["partition_coefficient"] > 0)
    ].copy()

    if valid.empty:
        ax.text(0.5, 0.5, "No valid data for Size vs K_part", ha="center", va="center", transform=ax.transAxes)
        return ax

    hue_col = "condensate_class" if "condensate_class" in valid.columns else None
    palette = {k: v for k, v in CLASS_COLOR_PALETTE.items() if k in valid[hue_col].values} if hue_col else None

    sns.scatterplot(
        data=valid,
        x=size_col,
        y="partition_coefficient",
        hue=hue_col,
        palette=palette,
        alpha=0.75,
        s=35,
        edgecolor="none",
        ax=ax,
    )

    #K = 1.0 baseline line
    ax.axhline(1.0, color="#d62728", linestyle="--", linewidth=1.2, label="K=1.0 (Nucleoplasm Baseline)")

    # Median K annotation
    med_k = valid["partition_coefficient"].median()
    ax.axhline(med_k, color="#2ca02c", linestyle=":", linewidth=1.2, label=f"Median K = {med_k:.2f}")

    ax.set_xscale("log")
    ax.set_xlabel(r"Condensate Volume ($\mu\mathrm{m}^3$)", fontsize=11, fontweight="bold")
    ax.set_ylabel(r"Partitioning Coefficient ($K_{\mathrm{part}}$)", fontsize=11, fontweight="bold")
    title = f"Size vs. Partitioning ($K_{{part}}$)"
    if sample_name:
        title += f" — {sample_name}"
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)

    if hue_col and ax.get_legend():
        ax.legend(title="Condensate Class", fontsize=8, title_fontsize=9, loc="upper right", framealpha=0.9)

    if standalone and output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(out_p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return ax


def plot_phase_diagram_anisotropy_vs_partitioning(
    df: pd.DataFrame,
    output_path: Path | str | None = None,
    ax: plt.Axes | None = None,
    sample_name: str | None = None,
) -> plt.Axes:
    """Plot 2D Biophysical Phase Diagram: Fractional Anisotropy (FA) vs. Partition Coefficient (K_part)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)

    fa_col = "A_object" if "A_object" in df.columns else None
    if not fa_col or "partition_coefficient" not in df.columns:
        if standalone:
            plt.close(fig)
        return ax

    valid = df[
        np.isfinite(df[fa_col])
        & np.isfinite(df["partition_coefficient"])
        & (df["partition_coefficient"] > 0)
    ].copy()

    if valid.empty:
        ax.text(0.5, 0.5, "No valid data for FA vs K_part", ha="center", va="center", transform=ax.transAxes)
        return ax

    max_k = max(10.0, float(valid["partition_coefficient"].quantile(0.99)) * 1.15)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.0, max_k)

    #Quadrant Shading
    #Top-Left: Unconstrained Droplet (Low FA, High K)
    ax.axvspan(0.0, 0.65, ymin=1.5 / max_k, ymax=1.0, color="#2ca02c", alpha=0.08)
    #Top-Right: Wetted / Chromatin-bound (High FA, High K)
    ax.axvspan(0.65, 1.0, ymin=1.5 / max_k, ymax=1.0, color="#1f77b4", alpha=0.08)
    #Bottom: Diffuse / Background (Low K)
    ax.axhspan(0.0, 1.5, color="#7f7f7f", alpha=0.08)

    #Quadrant Labels
    ax.text(0.32, max_k * 0.92, "Spherical LLPS Droplets\n(e.g., Free Speckles)",
            ha="center", va="center", fontsize=8.5, color="#2ca02c", fontweight="bold", alpha=0.85)
    ax.text(0.82, max_k * 0.92, "Wetted / Chromatin-Bound\n(e.g., POL II Hubs)",
            ha="center", va="center", fontsize=8.5, color="#1f77b4", fontweight="bold", alpha=0.85)
    ax.text(0.50, 0.75, "Diffuse / Transition Zone (K < 1.5)",
            ha="center", va="center", fontsize=8.5, color="#555555", alpha=0.8)

    #Threshold lines
    ax.axvline(0.65, color="#888888", linestyle=":", linewidth=1.0)
    ax.axhline(1.5, color="#888888", linestyle=":", linewidth=1.0)
    ax.axhline(1.0, color="#d62728", linestyle="--", linewidth=1.2, label="K=1.0 (Baseline)")

    hue_col = "condensate_class" if "condensate_class" in valid.columns else None
    palette = {k: v for k, v in CLASS_COLOR_PALETTE.items() if k in valid[hue_col].values} if hue_col else None

    sns.scatterplot(
        data=valid,
        x=fa_col,
        y="partition_coefficient",
        hue=hue_col,
        palette=palette,
        alpha=0.75,
        s=35,
        edgecolor="none",
        ax=ax,
    )

    ax.set_xlabel(r"3D Fractional Anisotropy ($FA_{\mathrm{object}}$)", fontsize=11, fontweight="bold")
    ax.set_ylabel(r"Partitioning Coefficient ($K_{\mathrm{part}}$)", fontsize=11, fontweight="bold")
    title = r"Biophysical Phase Diagram ($FA$ vs. $K_{\mathrm{part}}$)"
    if sample_name:
        title += f" — {sample_name}"
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)

    #Footnote about heuristic working thresholds
    ax.text(
        0.5, -0.15,
        "Note: Quadrant thresholds (K > 1.5, FA = 0.65 for 2:1 aspect ratio) represent recommended\nheuristic working boundaries for biological interpretation, not universal physical constants.",
        ha="center", va="top", transform=ax.transAxes, fontsize=7.5, color="#555555", style="italic"
    )

    if hue_col and ax.get_legend():
        ax.legend(title="Condensate Class", fontsize=8, title_fontsize=9, loc="upper right", framealpha=0.9)


    if standalone and output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(out_p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return ax


def export_partitioning_analysis(
    csv_or_df: Path | str | pd.DataFrame,
    output_dir: Path | str,
    file_stem: str | None = None,
) -> dict:
    """Generate both plots and summary metrics for Partitioning and Classification."""
    if isinstance(csv_or_df, (str, Path)):
        csv_path = Path(csv_or_df)
        df = pd.read_csv(csv_path)
        if file_stem is None:
            file_stem = csv_path.stem
    else:
        df = csv_or_df.copy()
        if file_stem is None:
            file_stem = "condensate_partitioning"

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    #Primary QC filtering if present
    plot_df = df[df["primary_qc_valid"]].copy() if "primary_qc_valid" in df.columns else df.copy()

    #Generate combined 2-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
    plot_size_vs_partitioning(plot_df, ax=ax1, sample_name=file_stem)
    plot_phase_diagram_anisotropy_vs_partitioning(plot_df, ax=ax2, sample_name=file_stem)
    plt.tight_layout()

    combined_png = out_dir / f"{file_stem}_partitioning_analysis.png"
    plt.savefig(combined_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    #Summary statistics
    k_vals = plot_df["partition_coefficient"].dropna() if "partition_coefficient" in plot_df.columns else pd.Series([], dtype=float)
    class_counts = plot_df["condensate_class"].value_counts(normalize=True) * 100 if "condensate_class" in plot_df.columns else pd.Series([], dtype=float)

    summary = {
        "sample_name": file_stem,
        "total_objects": len(plot_df),
        "median_partition_coefficient": float(k_vals.median()) if not k_vals.empty else np.nan,
        "mean_partition_coefficient": float(k_vals.mean()) if not k_vals.empty else np.nan,
        "class_percentages": class_counts.round(1).to_dict(),
        "plot_png": str(combined_png),
    }


    #Save summary JSON
    summary_json = out_dir / f"{file_stem}_partitioning_summary.json"
    with summary_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary
