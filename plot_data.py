import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import mplcursors
from pathlib import Path  

sns.set_theme(style="ticks", palette="muted")
plt.rcParams.update(
    {"font.sans-serif": "DejaVu Sans", "font.family": "sans-serif"}
)


csv_filename = "MED1_Output_Batch_3d.csv" 
csv_path = Path(csv_filename)

if not csv_path.exists():
    raise FileNotFoundError(f"File '{csv_filename}' not found in the current directory!")


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
    raise KeyError("In CSV there are no columns 'volume_bio_um3', 'area_bio_um2' or 'shape_metric_bio'!")

MAX_SIZE = 2.0  
MIN_SIZE = 0.0001 

df_filtered = df[(df["size_plot"] >= MIN_SIZE) & (df["size_plot"] <= MAX_SIZE)].copy()

if df_filtered.empty:
    raise ValueError(f"After filtering (MIN_SIZE={MIN_SIZE}, MAX_SIZE={MAX_SIZE}) No data left for plotting. Please check the input CSV file.")

counts_per_cell = (df_filtered.groupby(["filename", "cell_id"]).size().reset_index(name="condensate_count"))

fig, axes = plt.subplots(2, 2, figsize=(15, 11))

#Graph A: Distribution of Condensate Sizes
sns.histplot(
    df_filtered["size_plot"],
    kde=True,
    ax=axes[0, 0],
    color="#2b5c8f",
    bins=30 if len(df_filtered) > 30 else 10,
)
median_size = df_filtered["size_plot"].median()
axes[0, 0].axvline(
    median_size,
    color="#d9534f",
    linestyle="--",
    linewidth=2,
    label=f"Median: {median_size:.4f} {unit_label}",
)
axes[0, 0].set_title(title_a, fontsize=13, fontweight="bold")
axes[0, 0].set_xlabel(xlabel_text)
axes[0, 0].set_ylabel("Count")
axes[0, 0].legend()

#Graph B: Number of Condensates per Cell (FOV)
sns.boxplot(
    y=counts_per_cell["condensate_count"],
    ax=axes[0, 1],
    color="#8e44ad",
    width=0.3,
    boxprops=dict(alpha=0.7),
)
sns.stripplot(
    y=counts_per_cell["condensate_count"],
    ax=axes[0, 1],
    color="black",
    alpha=0.7,
    jitter=0.2,
    size=7,
)
axes[0, 1].set_title(
    "B) Number of Condensates per Cell (FOV)", fontsize=13, fontweight="bold"
)
axes[0, 1].set_ylabel("Number of Condensates")

fig, axes = plt.subplots(2, 3, figsize=(22, 12))

# Graph A: Distribution of Condensate Sizes
sns.histplot(
    df_filtered["size_plot"],
    kde=True,
    ax=axes[0, 0],
    color="#2b5c8f",
    bins=30 if len(df_filtered) > 30 else 10,
)
median_size = df_filtered["size_plot"].median()
axes[0, 0].axvline(
    median_size,
    color="#d9534f",
    linestyle="--",
    linewidth=2,
    label=f"Median: {median_size:.4f} {unit_label}",
)
axes[0, 0].set_title(title_a, fontsize=13, fontweight="bold")
axes[0, 0].set_xlabel(xlabel_text)
axes[0, 0].set_ylabel("Count")
axes[0, 0].legend()

# Graph B: Number of Condensates per Cell
sns.boxplot(
    y=counts_per_cell["condensate_count"],
    ax=axes[0, 1],
    color="#8e44ad",
    width=0.3,
    boxprops=dict(alpha=0.7),
)
sns.stripplot(
    y=counts_per_cell["condensate_count"],
    ax=axes[0, 1],
    color="black",
    alpha=0.7,
    jitter=0.2,
    size=7,
)
axes[0, 1].set_title(
    "B) Number of Condensates per Cell", fontsize=13, fontweight="bold"
)
axes[0, 1].set_ylabel("Number of Condensates")

# Graph C: Size vs. Integrated Density (Total Protein)
sns.scatterplot(
    data=df_filtered,
    x="size_plot",
    y="integrated_density",
    ax=axes[0, 2],
    color="#2e7d32",
    s=40,
    alpha=0.5,
)
sns.regplot(
    data=df_filtered,
    x="size_plot",
    y="integrated_density",
    ax=axes[0, 2],
    scatter=False,
    color="#1b5e20",
)
axes[0, 2].set_title(title_c, fontsize=13, fontweight="bold")
axes[0, 2].set_xlabel(xlabel_text)
axes[0, 2].set_ylabel("Integrated Density (a.u.)")

# Graph D: Signal Concentration (Violin Plot)
sns.violinplot(
    y=df_filtered["mean_intensity"],
    ax=axes[1, 0],
    color="#f57c00",
    inner="quartile",
    alpha=0.7
)
sns.stripplot(
    y=df_filtered["mean_intensity"],
    ax=axes[1, 0],
    color="black",
    alpha=0.2,
    jitter=0.15,
    size=3,
)
axes[1, 0].set_title(
    "D) Signal Concentration (Mean Intensity)", fontsize=13, fontweight="bold"
)
axes[1, 0].set_ylabel("Mean Intensity (a.u.)")

# Graph E: Size vs. Mean Intensity 
sns.scatterplot(
    data=df_filtered,
    x="size_plot",
    y="mean_intensity",
    ax=axes[1, 1],
    color="#c0392b",
    s=40,
    alpha=0.5,
)
sns.regplot(
    data=df_filtered,
    x="size_plot",
    y="mean_intensity",
    ax=axes[1, 1],
    scatter=False,
    color="#922b21",
)
axes[1, 1].set_title(
    "E) Size vs. Mean Intensity", fontsize=13, fontweight="bold"
)
axes[1, 1].set_xlabel(xlabel_text)
axes[1, 1].set_ylabel("Mean Intensity (a.u.)")

# Graph F: 2D Density (Cluster Topography)
sns.kdeplot(
    data=df_filtered,
    x="size_plot",
    y="mean_intensity",
    ax=axes[1, 2],
    fill=True,
    cmap="YlOrBr",
    thresh=0.05,
    alpha=0.8
)
axes[1, 2].set_title(
    "F) 2D Density (Cluster Topography)", fontsize=13, fontweight="bold"
)
axes[1, 2].set_xlabel(xlabel_text)
axes[1, 2].set_ylabel("Mean Intensity (a.u.)")

cursor = mplcursors.cursor([axes[0, 2].collections[0], axes[1, 1].collections[0]], hover=False) # hover=False znamená, že musíš kliknout. Můžeš změnit na True pro najetí myší.

@cursor.connect("add")
def on_add(sel):
    row = df_filtered.iloc[sel.index]
    sel.annotation.set_text(
        f"Photo: {row['filename']}\n"
        f"Cell ID: {row['cell_id']}\n"
        f"Volume: {row['size_plot']:.3f} µm³\n"
        f"Intensity: {row['mean_intensity']:.0f} a.u."
    )
    sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)

cursor_b = mplcursors.cursor(axes[0, 1].collections[-1], hover=False)

@cursor_b.connect("add")
def on_add_b(sel):
    row = counts_per_cell.iloc[sel.index]
    sel.annotation.set_text(
        f"Photo: {row['filename']}\n"
        f"Cell ID: {row['cell_id']}\n"
        f"Condensate Count: {row['condensate_count']}"
    )
    sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)

plt.tight_layout()

output_plot = f"{folder_name}_Analysis_Plots_6.png"
plt.savefig(output_plot, dpi=300)

print(f"Graphs were successfully created and saved as '{output_plot}'")
plt.show()
