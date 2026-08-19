import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def generate_excel_stats(csv_filename):

    csv_path = Path(csv_filename)

    if not csv_path.exists():
        print(f"Error: File '{csv_filename}' not found!")
        return

    folder_name = csv_path.resolve().parent.name
    output_excel = csv_path.parent / f"{folder_name}_Detailed_Stats.xlsx"
    output_csv = csv_path.parent / f"{folder_name}_All_Condensates_With_Diameters.csv"

    df = pd.read_csv(csv_path)

    if "volume_bio_um3" in df.columns:
        is_3d = True
        df["equivalent_diameter_um"] = ((6 * df["volume_bio_um3"]) / np.pi) ** (1 / 3)
        size_col = "volume_bio_um3"
    elif "area_bio_um2" in df.columns:
        is_3d = False
        df["equivalent_diameter_um"] = 2 * np.sqrt(df["area_bio_um2"] / np.pi)
        size_col = "area_bio_um2"
    elif "shape_metric_bio" in df.columns:
        is_3d = df["is_3d"].iloc[0] if "is_3d" in df.columns else True
        size_col = "shape_metric_bio"
        if is_3d:
            df["equivalent_diameter_um"] = ((6 * df["shape_metric_bio"]) / np.pi) ** (1 / 3)
        else:
            df["equivalent_diameter_um"] = 2 * np.sqrt(df["shape_metric_bio"] / np.pi)
    else:
        print("Error: Missing columns for size calculation!")
        return

    df["equivalent_diameter_nm"] = df["equivalent_diameter_um"] * 1000.0


    summary_stats = (
        df.groupby(["filename", "cell_id"])
        .agg(
            condensate_count=("object_id", "count"),
            mean_size=(size_col, "mean"),
            median_size=(size_col, "median"),
            mean_diameter_nm=("equivalent_diameter_nm", "mean"),
            median_diameter_nm=("equivalent_diameter_nm", "median"),
            mean_intensity=("mean_intensity", "mean"),
            integrated_density_sum=("integrated_density", "sum"),
        )
        .reset_index()
    )

    df.to_csv(output_csv, index=False)
    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        summary_stats.to_excel(writer, sheet_name="Per_Cell_Summary", index=False)
        df.to_excel(writer, sheet_name="All_Condensates", index=False)

    print(f"Excel statistics were generated: {output_excel}")
    print(f"CSV with all data was generated: {output_csv}")


def generate_plots(csv_filename, min_size=0.0001, max_size=2.0):
    sns.set_theme(style="ticks", palette="muted")
    plt.rcParams.update({"font.sans-serif": "DejaVu Sans", "font.family": "sans-serif"})

    csv_path = Path(csv_filename)
    if not csv_path.exists():
        print(f"Error: File '{csv_filename}' not found!")
        return

    folder_name = csv_path.resolve().parent.name
    df = pd.read_csv(csv_path)


    if "volume_bio_um3" in df.columns:
        df["size_plot"] = df["volume_bio_um3"]
        xlabel_text = "Biological Volume (µm³)"
        title_a = "A) Condensate Volume Distribution"
        title_c = "C) Volume vs. Total Intensity"
        unit_label = "µm³"
    elif "area_bio_um2" in df.columns:
        df["size_plot"] = df["area_bio_um2"]
        xlabel_text = "Biological Area (µm²)"
        title_a = "A) Condensate Area Distribution"
        title_c = "C) Area vs. Total Intensity"
        unit_label = "µm²"
    elif "shape_metric_bio" in df.columns:
        df["size_plot"] = df["shape_metric_bio"]
        is_3d_mode = df["is_3d"].iloc[0] if "is_3d" in df.columns else True
        unit_label = "µm³" if is_3d_mode else "µm²"
        dim_name = "Volume" if is_3d_mode else "Area"
        xlabel_text = f"Biological {dim_name} ({unit_label})"
        title_a = f"A) Condensate {dim_name} Distribution"
        title_c = f"C) {dim_name} vs. Total Intensity"
    else:
        print("Error: Missing columns for size calculation!")
        return

    df_filtered = df[(df["size_plot"] >= min_size) & (df["size_plot"] <= max_size)].copy()

    if df_filtered.empty:
        print("After filtration, no data left for plotting. Please check the input CSV file.")
        return

    counts_per_cell = df_filtered.groupby(["filename", "cell_id"]).size().reset_index(name="condensate_count")

    fig, axes = plt.subplots(2, 3, figsize=(22, 12))

    sns.histplot(df_filtered["size_plot"], kde=True, ax=axes[0, 0], color="#2b5c8f", bins=30 if len(df_filtered) > 30 else 10)
    median_size = df_filtered["size_plot"].median()
    axes[0, 0].axvline(median_size, color="#d9534f", linestyle="--", linewidth=2, label=f"Median: {median_size:.4f} {unit_label}")
    axes[0, 0].set_title(title_a, fontsize=13, fontweight="bold")
    axes[0, 0].set_xlabel(xlabel_text)
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].legend()

    sns.boxplot(y=counts_per_cell["condensate_count"], ax=axes[0, 1], color="#8e44ad", width=0.3, boxprops=dict(alpha=0.7))
    sns.stripplot(y=counts_per_cell["condensate_count"], ax=axes[0, 1], color="black", alpha=0.7, jitter=0.2, size=7)
    axes[0, 1].set_title("B) Number of Condensates per Cell", fontsize=13, fontweight="bold")
    axes[0, 1].set_ylabel("Number of Condensates")

    sns.scatterplot(data=df_filtered, x="size_plot", y="integrated_density", ax=axes[0, 2], color="#2e7d32", s=40, alpha=0.5)
    sns.regplot(data=df_filtered, x="size_plot", y="integrated_density", ax=axes[0, 2], scatter=False, color="#1b5e20")
    axes[0, 2].set_title(title_c, fontsize=13, fontweight="bold")
    axes[0, 2].set_xlabel(xlabel_text)
    axes[0, 2].set_ylabel("Integrated Density (a.u.)")

    sns.violinplot(y=df_filtered["mean_intensity"], ax=axes[1, 0], color="#f57c00", inner="quartile", alpha=0.7)
    sns.stripplot(y=df_filtered["mean_intensity"], ax=axes[1, 0], color="black", alpha=0.2, jitter=0.15, size=3)
    axes[1, 0].set_title("D) Signal Concentration (Mean Intensity)", fontsize=13, fontweight="bold")
    axes[1, 0].set_ylabel("Mean Intensity (a.u.)")

    sns.scatterplot(data=df_filtered, x="size_plot", y="mean_intensity", ax=axes[1, 1], color="#c0392b", s=40, alpha=0.5)
    sns.regplot(data=df_filtered, x="size_plot", y="mean_intensity", ax=axes[1, 1], scatter=False, color="#922b21")
    axes[1, 1].set_title("E) Size vs. Mean Intensity", fontsize=13, fontweight="bold")
    axes[1, 1].set_xlabel(xlabel_text)
    axes[1, 1].set_ylabel("Mean Intensity (a.u.)")

    sns.kdeplot(data=df_filtered, x="size_plot", y="mean_intensity", ax=axes[1, 2], fill=True, cmap="YlOrBr", thresh=0.05, alpha=0.8)
    axes[1, 2].set_title("F) 2D Density (Cluster Topography)", fontsize=13, fontweight="bold")
    axes[1, 2].set_xlabel(xlabel_text)
    axes[1, 2].set_ylabel("Mean Intensity (a.u.)")


    plt.tight_layout()
    output_plot = csv_path.parent / f"{folder_name}_Analysis_Plots.png"
    plt.savefig(output_plot, dpi=300)
    plt.close(fig)

    print(f"Graphs were successfuly generated to: '{output_plot}'")


def _as_bool_series(series):
    #Interpret CSV booleans safely after pandas reloads them as text or bool
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "on"})


def generate_rezim_a_plots(csv_filename):
    #Create a separate Core-Shell FA + QC report without changing standard ExQt plots
    sns.set_theme(style="whitegrid", palette="colorblind")
    plt.rcParams.update({"font.sans-serif": "DejaVu Sans", "font.family": "sans-serif"})

    csv_path = Path(csv_filename)
    if not csv_path.exists():
        print(f"Error: File '{csv_filename}' not found!")
        return None

    df = pd.read_csv(csv_path)
    required_columns = {
        "A_shell", "A_core", "Delta_A_core_shell", "mode_a_primary_include",
        "A_object_valid", "A_shell_valid", "A_middle_valid", "A_core_valid",
        "mode_a_qc_reason",
    }
    missing_columns = sorted(required_columns.difference(df.columns))
    if missing_columns:
        print("Rezim A plots skipped: missing columns: " + ", ".join(missing_columns))
        return None

    #Primary plots intentionally include only complete and QC-approved FA rows
    valid_flags = np.ones(len(df), dtype=bool)
    for column in ("A_object_valid", "A_shell_valid", "A_middle_valid", "A_core_valid"):
        valid_flags &= _as_bool_series(df[column]).to_numpy()

    numeric_columns = ["A_shell", "A_core", "Delta_A_core_shell"]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    complete_fa = valid_flags & df[numeric_columns].notna().all(axis=1).to_numpy()
    primary_include = _as_bool_series(df["mode_a_primary_include"]).to_numpy()
    primary = df.loc[complete_fa & primary_include].copy()

    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    fig.suptitle(
        "Rezim A: Core-Shell Fractional Anisotropy and QC Overview\n"
        "Core-Shell panels contain only complete, QC-approved objects.",
        fontsize=15,
        fontweight="bold",
    )

    def no_primary_data(ax, title):
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.text(
            0.5, 0.5,
            "No primary QC-valid objects\nfor this batch.",
            ha="center", va="center", transform=ax.transAxes, fontsize=12,
        )
        ax.set_xticks([])
        ax.set_yticks([])

    paired_ax = axes[0, 0]
    paired_ax.set_title("A) Paired Shell vs. Core FA", fontsize=12, fontweight="bold")
    if primary.empty:
        no_primary_data(paired_ax, "A) Paired Shell vs. Core FA")
    else:
        rng = np.random.default_rng(0)
        for _, row in primary.iterrows():
            paired_ax.plot([0, 1], [row["A_shell"], row["A_core"]], color="#6c757d", alpha=0.45, linewidth=1)
        paired_ax.scatter(rng.normal(0, 0.025, len(primary)), primary["A_shell"], color="#0f766e", label="Shell", zorder=3)
        paired_ax.scatter(rng.normal(1, 0.025, len(primary)), primary["A_core"], color="#c2410c", label="Core", zorder=3)
        paired_ax.set_xticks([0, 1], ["Shell", "Core"])
        paired_ax.set_ylabel("Fractional Anisotropy (FA)")
        paired_ax.set_ylim(0, 1)
        paired_ax.legend(frameon=True)

    scatter_ax = axes[0, 1]
    scatter_ax.set_title("B) Shell FA vs. Core FA", fontsize=12, fontweight="bold")
    if primary.empty:
        no_primary_data(scatter_ax, "B) Shell FA vs. Core FA")
    else:
        scatter_ax.scatter(primary["A_shell"], primary["A_core"], s=48, color="#2563eb", alpha=0.78)
        scatter_ax.plot([0, 1], [0, 1], linestyle="--", color="#6b7280", linewidth=1, label="Core = shell")
        scatter_ax.set_xlim(0, 1)
        scatter_ax.set_ylim(0, 1)
        scatter_ax.set_xlabel("Shell FA")
        scatter_ax.set_ylabel("Core FA")
        scatter_ax.legend(frameon=True)

    delta_ax = axes[1, 0]
    delta_ax.set_title("C) Delta A = Core FA − Shell FA", fontsize=12, fontweight="bold")
    if primary.empty:
        no_primary_data(delta_ax, "C) Delta A = Core FA − Shell FA")
    else:
        sns.histplot(primary["Delta_A_core_shell"], bins=min(20, max(5, len(primary))), kde=len(primary) >= 5,
                    color="#7c3aed", ax=delta_ax)
        delta_ax.axvline(0, color="#374151", linestyle="--", linewidth=1)
        delta_ax.set_xlabel("Delta A (geometric difference)")
        delta_ax.set_ylabel("Object count")

    qc_ax = axes[1, 1]
    qc_ax.set_title("D) Rezim A QC Funnel", fontsize=12, fontweight="bold")
    qc_counts = pd.Series(
        [len(df), int(complete_fa.sum()), len(primary)],
        index=["All segmented\nobjects", "Complete FA\nset", "Primary\nQC-valid"],
    )
    bars = qc_ax.bar(qc_counts.index, qc_counts.values, color=["#94a3b8", "#f59e0b", "#16a34a"])
    qc_ax.set_ylabel("Object count")
    for bar, value in zip(bars, qc_counts.values):
        qc_ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(value), ha="center", va="bottom")

    excluded_reasons = (
        df.loc[~(complete_fa & primary_include), "mode_a_qc_reason"]
        .fillna("")
        .str.split(";")
        .explode()
    )
    excluded_reasons = excluded_reasons[excluded_reasons != ""]
    if not excluded_reasons.empty:
        top_reasons = excluded_reasons.value_counts().head(3)
        reason_text = "Top QC reasons:\n" + "\n".join(f"{name}: {count}" for name, count in top_reasons.items())
        qc_ax.text(1.04, 0.96, reason_text, transform=qc_ax.transAxes, va="top", fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output_plot = csv_path.parent / f"{csv_path.parent.name}_Rezim_A_Plots.png"
    fig.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Rezim A plots were generated to: '{output_plot}'")
    return output_plot
