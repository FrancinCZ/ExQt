from __future__ import annotations

import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns


CLASS_COLOR_PALETTE = {
    "Globular (Low FA) / Core-Enriched": "#2ca02c",
    "Globular (Low FA) / Shell-Enriched": "#ff7f0e",
    "Elongated (High FA) / Core-Enriched": "#1f77b4",
    "Elongated (High FA) / Shell-Enriched": "#d62728",
    "Punctate / Small Cluster": "#7f7f7f",
    "Unclassified": "#9467bd",
    #Backward compatibility aliases
    "Unconstrained Droplet (Spherical LLPS)": "#2ca02c",
    "Wetted / Chromatin-Bound Condensate": "#1f77b4",
    "Hollow / Shell-Dominated": "#ff7f0e",
    "Complex / Multi-Phase Aggregate": "#d62728",
}


def _get_particle_class_col(df: pd.DataFrame) -> str | None:
    if "particle_class" in df.columns:
        return "particle_class"
    if "condensate_class" in df.columns:
        return "condensate_class"
    return None


def plot_size_vs_partitioning(
    df: pd.DataFrame,
    output_path: Path | str | None = None,
    ax: plt.Axes | None = None,
    sample_name: str | None = None,
    show_legend: bool = True,
) -> plt.Axes:
    """Plot Particle Size (Volume in um3) vs. Partition Coefficient (K_part)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)

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

    hue_col = _get_particle_class_col(valid)
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
        legend=show_legend,
    )

    #K = 1.0 baseline line
    ax.axhline(1.0, color="#d62728", linestyle="--", linewidth=1.2, label="K=1.0 (Baseline)")

    #Median K annotation
    med_k = valid["partition_coefficient"].median()
    ax.axhline(med_k, color="#2ca02c", linestyle=":", linewidth=1.2, label=f"Median K = {med_k:.2f}")

    ax.set_xscale("log")
    ax.set_xlabel(r"Particle Volume ($\mu\mathrm{m}^3$)", fontsize=10, fontweight="bold")
    ax.set_ylabel(r"Partitioning ($K_{\mathrm{part}}$)", fontsize=10, fontweight="bold")
    title = r"Size vs. Partitioning ($K_{\mathrm{part}}$)"
    if standalone and sample_name:
        title += f" — {sample_name}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)

    if show_legend and hue_col and ax.get_legend():
        ax.legend(title="Class", fontsize=7.5, title_fontsize=8.5, loc="upper right", framealpha=0.9)
    elif not show_legend and ax.get_legend():
        ax.get_legend().remove()

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
    show_legend: bool = True,
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
    ax.axvspan(0.0, 0.65, ymin=1.5 / max_k, ymax=1.0, color="#2ca02c", alpha=0.08)
    ax.axvspan(0.65, 1.0, ymin=1.5 / max_k, ymax=1.0, color="#1f77b4", alpha=0.08)
    ax.axhspan(0.0, 1.5, color="#7f7f7f", alpha=0.08)

    #Objective Quadrant Labels
    ax.text(0.32, max_k * 0.94, "Globular (Low FA)",
            ha="center", va="center", fontsize=8.5, color="#2ca02c", fontweight="bold", alpha=0.85)
    ax.text(0.82, max_k * 0.94, "Elongated (High FA)",
            ha="center", va="center", fontsize=8.5, color="#1f77b4", fontweight="bold", alpha=0.85)
    ax.text(0.50, 0.75, "Low Contrast (K < 1.5)",
            ha="center", va="center", fontsize=8.5, color="#555555", alpha=0.8)

    #Threshold lines
    ax.axvline(0.65, color="#888888", linestyle=":", linewidth=1.0, label="Globular Boundary (0.65)")
    ax.axhline(1.5, color="#888888", linestyle=":", linewidth=1.0)
    ax.axhline(1.0, color="#d62728", linestyle="--", linewidth=1.2, label="K=1.0 (Baseline)")

    hue_col = _get_particle_class_col(valid)
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
        legend=show_legend,
    )

    ax.set_xlabel(r"3D Fractional Anisotropy ($FA_{\mathrm{object}}$)", fontsize=10, fontweight="bold")
    ax.set_ylabel(r"Partitioning ($K_{\mathrm{part}}$)", fontsize=10, fontweight="bold")
    title = r"Phase Diagram ($FA$ vs. $K_{\mathrm{part}}$)"
    if standalone and sample_name:
        title += f" — {sample_name}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)

    if show_legend and hue_col and ax.get_legend():
        ax.legend(title="Class", fontsize=7.5, title_fontsize=8.5, loc="upper right", framealpha=0.9)
    elif not show_legend and ax.get_legend():
        ax.get_legend().remove()

    if standalone and output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(out_p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return ax


def plot_volume_vs_intensity(
    df: pd.DataFrame,
    output_path: Path | str | None = None,
    ax: plt.Axes | None = None,
    sample_name: str | None = None,
    show_legend: bool = True,
) -> plt.Axes:
    """Plot Particle Volume vs. Mean Internal Intensity (Test of Phase Coexistence / Lever Rule)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)

    size_col = "volume_bio_um3" if "volume_bio_um3" in df.columns else "shape_metric_bio"
    int_col = "mean_intensity" if "mean_intensity" in df.columns else "intensity_mean"
    if size_col not in df.columns or int_col not in df.columns:
        if standalone:
            plt.close(fig)
        return ax

    valid = df[
        np.isfinite(df[size_col])
        & np.isfinite(df[int_col])
        & (df[size_col] > 0)
        & (df[int_col] > 0)
    ].copy()

    if valid.empty:
        ax.text(0.5, 0.5, "No valid data for Volume vs Intensity", ha="center", va="center", transform=ax.transAxes)
        return ax

    hue_col = _get_particle_class_col(valid)
    palette = {k: v for k, v in CLASS_COLOR_PALETTE.items() if k in valid[hue_col].values} if hue_col else None

    sns.scatterplot(
        data=valid,
        x=size_col,
        y=int_col,
        hue=hue_col,
        palette=palette,
        alpha=0.75,
        s=35,
        edgecolor="none",
        ax=ax,
        legend=show_legend,
    )

    #Median internal intensity reference line
    med_int = valid[int_col].median()
    ax.axhline(med_int, color="#333333", linestyle="--", linewidth=1.2, label=f"Median Density ({med_int:.1f} ADU)")

    ax.set_xscale("log")
    ax.set_xlabel(r"Particle Volume ($\mu\mathrm{m}^3$)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Mean Intensity (ADU)", fontsize=10, fontweight="bold")
    title = "Volume vs. Internal Density (Phase Coexistence)"
    if standalone and sample_name:
        title += f" — {sample_name}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)

    if show_legend and hue_col and ax.get_legend():
        ax.legend(title="Class", fontsize=7.5, title_fontsize=8.5, loc="upper right", framealpha=0.9)
    elif not show_legend and ax.get_legend():
        ax.get_legend().remove()

    if standalone and output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(out_p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return ax


def plot_volume_vs_anisotropy(
    df: pd.DataFrame,
    output_path: Path | str | None = None,
    ax: plt.Axes | None = None,
    sample_name: str | None = None,
    show_legend: bool = True,
) -> plt.Axes:
    """Plot Particle Volume vs. 3D Fractional Anisotropy (Test of Surface Tension vs. Optical Limit)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)

    size_col = "volume_bio_um3" if "volume_bio_um3" in df.columns else "shape_metric_bio"
    fa_col = "A_object" if "A_object" in df.columns else None
    if size_col not in df.columns or not fa_col:
        if standalone:
            plt.close(fig)
        return ax

    valid = df[
        np.isfinite(df[size_col])
        & np.isfinite(df[fa_col])
        & (df[size_col] > 0)
    ].copy()

    if valid.empty:
        ax.text(0.5, 0.5, "No valid data for Volume vs FA", ha="center", va="center", transform=ax.transAxes)
        return ax

    hue_col = _get_particle_class_col(valid)
    palette = {k: v for k, v in CLASS_COLOR_PALETTE.items() if k in valid[hue_col].values} if hue_col else None

    sns.scatterplot(
        data=valid,
        x=size_col,
        y=fa_col,
        hue=hue_col,
        palette=palette,
        alpha=0.75,
        s=35,
        edgecolor="none",
        ax=ax,
        legend=show_legend,
    )

    #Reference lines for PSF limit and Globular boundary
    ax.axhline(0.82, color="#d62728", linestyle=":", linewidth=1.2, label="PSF Elongation Limit (~0.82)")
    ax.axhline(0.65, color="#888888", linestyle="--", linewidth=1.2, label="Globular Boundary (0.65)")

    ax.set_xscale("log")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"Particle Volume ($\mu\mathrm{m}^3$)", fontsize=10, fontweight="bold")
    ax.set_ylabel(r"3D Fractional Anisotropy ($FA_{\mathrm{object}}$)", fontsize=10, fontweight="bold")
    title = "Volume vs. Shape Anisotropy (Surface Relaxation)"
    if standalone and sample_name:
        title += f" — {sample_name}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)

    if show_legend and hue_col and ax.get_legend():
        ax.legend(title="Class", fontsize=7.5, title_fontsize=8.5, loc="upper right", framealpha=0.9)
    elif not show_legend and ax.get_legend():
        ax.get_legend().remove()

    if standalone and output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(out_p, dpi=300, bbox_inches="tight")
        plt.close(fig)

    return ax


def export_unified_pdf_report(
    df: pd.DataFrame,
    pdf_path: Path | str,
    sample_name: str = "Sample",
    min_size: float | None = None,
    max_size: float | None = None,
) -> None:
    """Generate a clean, publication-ready 2-page vector PDF report summarizing all biophysical metrics."""
    from matplotlib.lines import Line2D

    out_pdf = Path(pdf_path)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    plot_df = df[df["primary_qc_valid"]].copy() if "primary_qc_valid" in df.columns else df.copy()

    #Apply size filter consistent with GUI preview
    size_col = "volume_bio_um3" if "volume_bio_um3" in plot_df.columns else ("area_bio_um2" if "area_bio_um2" in plot_df.columns else None)
    if size_col is not None:
        if min_size is not None:
            plot_df = plot_df[plot_df[size_col] >= min_size]
        if max_size is not None:
            plot_df = plot_df[plot_df[size_col] <= max_size]

    with PdfPages(out_pdf) as pdf:
        fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11.5), dpi=300)
        
        plot_phase_diagram_anisotropy_vs_partitioning(plot_df, ax=ax1, show_legend=False)
        plot_size_vs_partitioning(plot_df, ax=ax2, show_legend=False)
        plot_volume_vs_intensity(plot_df, ax=ax3, show_legend=False)
        plot_volume_vs_anisotropy(plot_df, ax=ax4, show_legend=False)

        # 1. Class handles (colored dots)
        class_handles, class_labels = [], []
        class_col = _get_particle_class_col(plot_df)
        if class_col and class_col in plot_df.columns:
            present = [c for c in CLASS_COLOR_PALETTE if c in plot_df[class_col].values]
            for c in present:
                class_handles.append(
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=CLASS_COLOR_PALETTE[c], markersize=8, markeredgecolor="none")
                )
                class_labels.append(c)

        #2. Reference lines
        k_vals = plot_df["partition_coefficient"].dropna() if "partition_coefficient" in plot_df.columns else pd.Series([], dtype=float)
        int_col = "mean_intensity" if "mean_intensity" in plot_df.columns else "intensity_mean"
        med_int_val = plot_df[int_col].median() if int_col in plot_df.columns and not plot_df[int_col].dropna().empty else None

        ref_handles = [
            Line2D([0], [0], color="#d62728", linestyle="--", linewidth=1.4),
            Line2D([0], [0], color="#2ca02c", linestyle=":", linewidth=1.4),
            Line2D([0], [0], color="#333333", linestyle="--", linewidth=1.4),
            Line2D([0], [0], color="#d62728", linestyle=":", linewidth=1.4),
            Line2D([0], [0], color="#888888", linestyle="--", linewidth=1.4),
        ]
        ref_labels = [
            "K = 1.0 (Baseline)",
            f"Median K = {k_vals.median():.2f}" if not k_vals.empty else "Median K",
            f"Median Density = {med_int_val:.0f} ADU" if med_int_val is not None else "Median Density",
            "PSF Limit (~0.82)",
            "Globular Boundary (0.65)",
        ]

        all_handles = class_handles + ref_handles
        all_labels = class_labels + ref_labels

        if all_handles:
            fig1.legend(
                all_handles, all_labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.015),
                ncol=min(len(all_labels), 4),
                fontsize=8.5,
                frameon=True,
                facecolor="#fcfcfc",
                edgecolor="#cccccc",
                framealpha=0.98,
                title="Particle Classes & Biophysical Reference Lines",
                title_fontsize=9,
            )

        fig1.suptitle(f"ExQt Biophysical Particle Profiling — {sample_name}", fontsize=14, fontweight="bold", y=0.98)
        plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.17, wspace=0.25, hspace=0.32)
        pdf.savefig(fig1)
        plt.close(fig1)

        fig2, (ax_hist, ax_table) = plt.subplots(2, 1, figsize=(14, 11), dpi=300, gridspec_kw={"height_ratios": [1.1, 1]})

        #Panel 1: Size Distribution Histogram
        size_col = "volume_bio_um3" if "volume_bio_um3" in plot_df.columns else "shape_metric_bio"
        if size_col in plot_df.columns:
            sizes = plot_df[size_col].dropna()
            sizes = sizes[sizes > 0]
            if not sizes.empty:
                sns.histplot(sizes, log_scale=True, kde=True, color="#1f77b4", ax=ax_hist, bins=30)
                ax_hist.set_title("Particle Volume Distribution", fontsize=11, fontweight="bold")
                ax_hist.set_xlabel(r"Biological Volume ($\mu\mathrm{m}^3$)", fontsize=10, fontweight="bold")
                ax_hist.set_ylabel("Count", fontsize=10, fontweight="bold")
                ax_hist.grid(True, linestyle="--", alpha=0.4)

        #Panel 2: Summary Metrics Table
        ax_table.axis("off")
        class_counts = plot_df[class_col].value_counts() if class_col else pd.Series([], dtype=int)

        table_data = [
            ["Metric", "Value"],
            ["Sample ID", str(sample_name)],
            ["Total Primary Particles", f"{len(plot_df):,}"],
            ["Median Partitioning (K_part)", f"{k_vals.median():.2f}" if not k_vals.empty else "N/A"],
            ["Mean Partitioning (K_part)", f"{k_vals.mean():.2f}" if not k_vals.empty else "N/A"],
            ["Median Particle Volume", f"{plot_df[size_col].median():.4f} µm³" if size_col in plot_df.columns else "N/A"],
        ]
        for cname, count in class_counts.items():
            pct = (count / len(plot_df)) * 100 if len(plot_df) > 0 else 0
            table_data.append([f"Class: {cname}", f"{count:,} ({pct:.1f}%)"])

        table = ax_table.table(
            cellText=table_data,
            loc="center",
            cellLoc="left",
            colWidths=[0.5, 0.4],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9.5)
        table.scale(1.0, 1.4)
        for i in range(2):
            table[(0, i)].set_facecolor("#2c3e50")
            table[(0, i)].get_text().set_color("white")
            table[(0, i)].get_text().set_fontweight("bold")

        fig2.suptitle(f"Morphological Distribution & Summary — {sample_name}", fontsize=14, fontweight="bold", y=0.97)
        plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.06, hspace=0.32)
        pdf.savefig(fig2)
        plt.close(fig2)


def export_partitioning_analysis(
    csv_or_df: Path | str | pd.DataFrame,
    output_dir: Path | str,
    file_stem: str | None = None,
    export_pdf: bool = True,
    export_pngs: bool = True,
    min_size: float | None = None,
    max_size: float | None = None,
) -> dict:
    """Generate plots, executive 2-page PDF report, and summary metrics."""
    from matplotlib.lines import Line2D

    if isinstance(csv_or_df, (str, Path)):
        csv_path = Path(csv_or_df)
        df = pd.read_csv(csv_path)
        if file_stem is None:
            file_stem = csv_path.stem
    else:
        df = csv_or_df.copy()
        if file_stem is None:
            file_stem = "particle_analysis"

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_df = df[df["primary_qc_valid"]].copy() if "primary_qc_valid" in df.columns else df.copy()

    #Apply size filter consistent with GUI preview (volume_bio_um3 for 3D, area_bio_um2 for 2D)
    size_col = "volume_bio_um3" if "volume_bio_um3" in plot_df.columns else ("area_bio_um2" if "area_bio_um2" in plot_df.columns else None)
    if size_col is not None:
        if min_size is not None:
            plot_df = plot_df[plot_df[size_col] >= min_size]
        if max_size is not None:
            plot_df = plot_df[plot_df[size_col] <= max_size]

    #Compute class column and K series early so summary statistics work regardless of export_pngs
    class_col = _get_particle_class_col(plot_df)
    k_vals = plot_df["partition_coefficient"].dropna() if "partition_coefficient" in plot_df.columns else pd.Series([], dtype=float)

    #Generate 4-panel combined figure (PNG) with unified bottom legend and clean spacing
    combined_png = out_dir / f"{file_stem}_partitioning_analysis.png"
    if export_pngs:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11.5), dpi=300)
        plot_phase_diagram_anisotropy_vs_partitioning(plot_df, ax=ax1, show_legend=False)
        plot_size_vs_partitioning(plot_df, ax=ax2, show_legend=False)
        plot_volume_vs_intensity(plot_df, ax=ax3, show_legend=False)
        plot_volume_vs_anisotropy(plot_df, ax=ax4, show_legend=False)

        # 1. Class handles (colored dots)
        class_handles, class_labels = [], []
        if class_col and class_col in plot_df.columns:
            present = [c for c in CLASS_COLOR_PALETTE if c in plot_df[class_col].values]
            for c in present:
                class_handles.append(
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=CLASS_COLOR_PALETTE[c], markersize=8, markeredgecolor="none")
                )
                class_labels.append(c)

        #2. Reference lines
        int_col = "mean_intensity" if "mean_intensity" in plot_df.columns else "intensity_mean"
        med_int_val = plot_df[int_col].median() if int_col in plot_df.columns and not plot_df[int_col].dropna().empty else None

        ref_handles = [
            Line2D([0], [0], color="#d62728", linestyle="--", linewidth=1.4),
            Line2D([0], [0], color="#2ca02c", linestyle=":", linewidth=1.4),
            Line2D([0], [0], color="#333333", linestyle="--", linewidth=1.4),
            Line2D([0], [0], color="#d62728", linestyle=":", linewidth=1.4),
            Line2D([0], [0], color="#888888", linestyle="--", linewidth=1.4),
        ]
        ref_labels = [
            "K = 1.0 (Baseline)",
            f"Median K = {k_vals.median():.2f}" if not k_vals.empty else "Median K",
            f"Median Density = {med_int_val:.0f} ADU" if med_int_val is not None else "Median Density",
            "PSF Limit (~0.82)",
            "Globular Boundary (0.65)",
        ]

        all_handles = class_handles + ref_handles
        all_labels = class_labels + ref_labels

        if all_handles:
            fig.legend(
                all_handles, all_labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.015),
                ncol=min(len(all_labels), 4),
                fontsize=8.5,
                frameon=True,
                facecolor="#fcfcfc",
                edgecolor="#cccccc",
                framealpha=0.98,
                title="Particle Classes & Biophysical Reference Lines",
                title_fontsize=9,
            )

        fig.suptitle(f"ExQt Biophysical Particle Profiling — {file_stem}", fontsize=14, fontweight="bold", y=0.98)
        plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.17, wspace=0.25, hspace=0.32)
        plt.savefig(combined_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

    #Generate unified 2-page PDF report
    pdf_report_path = out_dir / f"{file_stem}_ExQt_Report.pdf"
    if export_pdf:
        export_unified_pdf_report(plot_df, pdf_report_path, sample_name=file_stem)

    #Summary statistics
    class_counts = plot_df[class_col].value_counts(normalize=True) * 100 if class_col else pd.Series([], dtype=float)

    summary = {
        "sample_name": file_stem,
        "total_particles": len(plot_df),
        "total_objects": len(plot_df),  # backward compatibility alias
        "median_partition_coefficient": float(k_vals.median()) if not k_vals.empty else np.nan,
        "mean_partition_coefficient": float(k_vals.mean()) if not k_vals.empty else np.nan,
        "class_percentages": class_counts.round(1).to_dict(),
        "plot_png": str(combined_png) if export_pngs else None,
        "report_pdf": str(pdf_report_path) if export_pdf else None,
    }

    #Save summary JSON
    summary_json = out_dir / f"{file_stem}_partitioning_summary.json"
    with summary_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary
