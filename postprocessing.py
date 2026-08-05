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
        df.groupby("filename")
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

    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        df.to_csv(output_csv, index=False)
        summary_stats.to_excel(writer, sheet_name="Per_Cell_Summary", index=False)
        df.to_excel(writer, sheet_name="All_Condensates", index=False)

    print(f"Excel statistics were generated: {output_excel}")
    print(f"CSV with all data was generated: {output_csv}")


def generate_plots(csv_filename):
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

    MAX_SIZE = 2.0
    MIN_SIZE = 0.0001
    df_filtered = df[(df["size_plot"] >= MIN_SIZE) & (df["size_plot"] <= MAX_SIZE)].copy()

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


if __name__ == "__main__":
    generate_plots("Magnify_test_Output_Batch_3d.csv")