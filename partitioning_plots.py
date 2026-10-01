from __future__ import annotations

import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns

from reference_values import reference_value


#Legend text for the reference lines; the values are working values (reference_values.py).
def _reference_labels():
    return {
        "psf": f"PSF-like elongation, FA {reference_value('fa_psf_elongation_limit'):g} (working value)",
        "globular": f"Globular/Elongated boundary, FA {reference_value('fa_globular_boundary'):g} (working value)",
        "low_k": f"K < {reference_value('k_low_contrast'):g} (working value)",
    }


CLASS_COLOR_PALETTE = {
    "Globular (Low FA) / Core-Enriched": "#2ca02c",
    "Globular (Low FA) / Shell-Enriched": "#ff7f0e",
    "Elongated (High FA) / Core-Enriched": "#1f77b4",
    "Elongated (High FA) / Shell-Enriched": "#d62728",
    "Globular (Low FA) / No Significant Gradient": "#98df8a",
    "Elongated (High FA) / No Significant Gradient": "#aec7e8",
    "Globular (Low FA) / Profile-Unavailable": "#c7c7c7",
    "Elongated (High FA) / Profile-Unavailable": "#c7c7c7",
    "Punctate / Small Cluster": "#7f7f7f",
    "Unclassified": "#9467bd",
    #Backward compatibility aliases
    "Unconstrained Droplet (Spherical LLPS)": "#2ca02c",
    "Wetted / Chromatin-Bound Condensate": "#1f77b4",
    "Hollow / Shell-Dominated": "#ff7f0e",
    "Complex / Multi-Phase Aggregate": "#d62728",
}


PSF_INTERPRETATION_WARNING = (
    "Core-Enriched is the expected result of PSF blurring even for a homogeneous object; "
    "it is not evidence of biological core enrichment until compared with a null/PSF simulation. "
    "Classes use a two-sided Welch t-test core vs shell voxels (alpha = 0.01); only QC-passed "
    "(primary) objects are shown and K_part only where K_valid."
)


#Return the QC-passed objects for reporting, with K masked where K_valid is False.
def _report_subset(df: pd.DataFrame) -> pd.DataFrame:
    if "primary_qc_valid" not in df.columns:
        raise ValueError(
            "Partitioning report needs the QC-annotated frame with 'primary_qc_valid' "
            "(postprocessing._prepare_reporting_frames), not the raw batch CSV."
        )
    passed = df["primary_qc_valid"].astype(str).str.strip().str.lower().isin({"true", "1"})
    subset = df[passed].copy()
    if "partition_coefficient" in subset.columns and "K_valid" in subset.columns:
        k_valid = subset["K_valid"].astype(str).str.strip().str.lower().isin({"true", "1"})
        subset.loc[~k_valid, "partition_coefficient"] = np.nan
    return subset


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

    ax.axhline(1.0, color="#d62728", linestyle="--", linewidth=1.2, label="K=1.0 (Baseline)")

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

    fa_boundary = reference_value("fa_globular_boundary")
    k_low = reference_value("k_low_contrast")
    labels = _reference_labels()

    #Quadrant shading from the central working values
    ax.axvspan(0.0, fa_boundary, ymin=k_low / max_k, ymax=1.0, color="#2ca02c", alpha=0.08)
    ax.axvspan(fa_boundary, 1.0, ymin=k_low / max_k, ymax=1.0, color="#1f77b4", alpha=0.08)
    ax.axhspan(0.0, k_low, color="#7f7f7f", alpha=0.08)

    ax.text(fa_boundary / 2, max_k * 0.94, "Globular (Low FA)",
            ha="center", va="center", fontsize=8.5, color="#2ca02c", fontweight="bold", alpha=0.85)
    ax.text((fa_boundary + 1.0) / 2, max_k * 0.94, "Elongated (High FA)",
            ha="center", va="center", fontsize=8.5, color="#1f77b4", fontweight="bold", alpha=0.85)
    ax.text(0.50, k_low / 2, labels["low_k"],
            ha="center", va="center", fontsize=8.5, color="#555555", alpha=0.8)

    ax.axvline(fa_boundary, color="#888888", linestyle=":", linewidth=1.0, label=labels["globular"])
    ax.axhline(k_low, color="#888888", linestyle=":", linewidth=1.0)
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
    """Plot Particle Volume vs. Mean Internal Intensity (descriptive scatter, not a test of any model)."""
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

    med_int = valid[int_col].median()
    ax.axhline(med_int, color="#333333", linestyle="--", linewidth=1.2, label=f"Median Density ({med_int:.1f} ADU)")

    ax.set_xscale("log")
    ax.set_xlabel(r"Particle Volume ($\mu\mathrm{m}^3$)", fontsize=10, fontweight="bold")
    ax.set_ylabel("Mean Intensity (ADU)", fontsize=10, fontweight="bold")
    title = "Volume vs. Mean Intensity"
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

    #Reference lines (working values from reference_values.py)
    labels = _reference_labels()
    ax.axhline(reference_value("fa_psf_elongation_limit"), color="#d62728", linestyle=":", linewidth=1.2, label=labels["psf"])
    ax.axhline(reference_value("fa_globular_boundary"), color="#888888", linestyle="--", linewidth=1.2, label=labels["globular"])

    ax.set_xscale("log")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"Particle Volume ($\mu\mathrm{m}^3$)", fontsize=10, fontweight="bold")
    ax.set_ylabel(r"3D Fractional Anisotropy ($FA_{\mathrm{object}}$)", fontsize=10, fontweight="bold")
    title = "Volume vs. Shape Anisotropy"
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


def _size_filtered(df, min_size=None, max_size=None):
    """QC-passed objects within the GUI size range (volume in 3D, area in 2D)."""
    plot_df = _report_subset(df)
    size_col = next((c for c in ("volume_bio_um3", "area_bio_um2") if c in plot_df.columns), None)
    if size_col is not None:
        if min_size is not None:
            plot_df = plot_df[plot_df[size_col] >= min_size]
        if max_size is not None:
            plot_df = plot_df[plot_df[size_col] <= max_size]
    return plot_df


def _k_values(plot_df):
    return plot_df["partition_coefficient"].dropna() if "partition_coefficient" in plot_df.columns else pd.Series([], dtype=float)


def _four_panel_figure(plot_df, title):
    """The four summary plots with one shared legend (classes + reference lines)."""
    from matplotlib.lines import Line2D

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11.5), dpi=300)
    plot_phase_diagram_anisotropy_vs_partitioning(plot_df, ax=ax1, show_legend=False)
    plot_size_vs_partitioning(plot_df, ax=ax2, show_legend=False)
    plot_volume_vs_intensity(plot_df, ax=ax3, show_legend=False)
    plot_volume_vs_anisotropy(plot_df, ax=ax4, show_legend=False)

    class_col = _get_particle_class_col(plot_df)
    classes = [c for c in CLASS_COLOR_PALETTE if class_col and c in plot_df[class_col].values]
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CLASS_COLOR_PALETTE[c], markersize=8,
                    markeredgecolor="none") for c in classes]
    k_vals = _k_values(plot_df)
    int_col = "mean_intensity" if "mean_intensity" in plot_df.columns else "intensity_mean"
    med_int = plot_df[int_col].median() if int_col in plot_df.columns and not plot_df[int_col].dropna().empty else None
    handles += [Line2D([0], [0], color=color, linestyle=style, linewidth=1.4) for color, style in (
        ("#d62728", "--"), ("#2ca02c", ":"), ("#333333", "--"), ("#d62728", ":"), ("#888888", "--"))]
    labels = classes + [
        "K = 1.0 (Baseline)",
        f"Median K = {k_vals.median():.2f}" if not k_vals.empty else "Median K",
        f"Median Density = {med_int:.0f} ADU" if med_int is not None else "Median Density",
        _reference_labels()["psf"],
        _reference_labels()["globular"],
    ]
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.015), ncol=min(len(labels), 4),
            fontsize=8.5, frameon=True, facecolor="#fcfcfc", edgecolor="#cccccc", framealpha=0.98,
            title="Particle Classes & Biophysical Reference Lines", title_fontsize=9)
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)
    plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.17, wspace=0.25, hspace=0.32)
    return fig


def export_unified_pdf_report(
    df: pd.DataFrame,
    pdf_path: Path | str,
    sample_name: str = "Sample",
    min_size: float | None = None,
    max_size: float | None = None,
) -> None:
    """2-page PDF: the four summary plots, then size distribution and a summary table."""
    out_pdf = Path(pdf_path)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plot_df = _size_filtered(df, min_size, max_size)

    with PdfPages(out_pdf) as pdf:
        fig1 = _four_panel_figure(plot_df, f"ExQt Biophysical Particle Profiling — {sample_name}")
        pdf.savefig(fig1)
        plt.close(fig1)

        fig2, (ax_hist, ax_table) = plt.subplots(2, 1, figsize=(14, 11), dpi=300, gridspec_kw={"height_ratios": [1.1, 1]})
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

        ax_table.axis("off")
        class_col = _get_particle_class_col(plot_df)
        class_counts = plot_df[class_col].value_counts() if class_col else pd.Series([], dtype=int)
        k_vals = _k_values(plot_df)
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
        table = ax_table.table(cellText=table_data, loc="center", cellLoc="left", colWidths=[0.5, 0.4])
        table.auto_set_font_size(False)
        table.set_fontsize(9.5)
        table.scale(1.0, 1.4)
        for i in range(2):
            table[(0, i)].set_facecolor("#2c3e50")
            table[(0, i)].get_text().set_color("white")
            table[(0, i)].get_text().set_fontweight("bold")

        fig2.text(0.5, 0.015, PSF_INTERPRETATION_WARNING, ha="center", va="bottom", fontsize=8.5,
                style="italic", wrap=True, color="#8b0000")
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
    """Four-panel PNG, 2-page PDF and a summary JSON for the QC-passed objects."""
    if isinstance(csv_or_df, (str, Path)):
        df = pd.read_csv(csv_or_df)
        file_stem = file_stem or Path(csv_or_df).stem
    else:
        df = csv_or_df.copy()
        file_stem = file_stem or "particle_analysis"

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_df = _size_filtered(df, min_size, max_size)
    class_col = _get_particle_class_col(plot_df)
    k_vals = _k_values(plot_df)

    combined_png = out_dir / f"{file_stem}_partitioning_analysis.png"
    if export_pngs:
        fig = _four_panel_figure(plot_df, f"ExQt Biophysical Particle Profiling — {file_stem}")
        plt.savefig(combined_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

    pdf_report_path = out_dir / f"{file_stem}_ExQt_Report.pdf"
    if export_pdf:
        export_unified_pdf_report(plot_df, pdf_report_path, sample_name=file_stem)

    class_counts = plot_df[class_col].value_counts(normalize=True) * 100 if class_col else pd.Series([], dtype=float)
    summary = {
        "sample_name": file_stem,
        "total_particles": len(plot_df),
        "total_objects": len(plot_df),  # backward compatibility alias
        "median_partition_coefficient": float(k_vals.median()) if not k_vals.empty else np.nan,
        "mean_partition_coefficient": float(k_vals.mean()) if not k_vals.empty else np.nan,
        "class_percentages": class_counts.round(1).to_dict(),
        "interpretation_warning": PSF_INTERPRETATION_WARNING,
        "plot_png": str(combined_png) if export_pngs else None,
        "report_pdf": str(pdf_report_path) if export_pdf else None,
    }
    with (out_dir / f"{file_stem}_partitioning_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    return summary
