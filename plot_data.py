import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path  

sns.set_theme(style="ticks", palette="muted")
plt.rcParams.update(
    {"font.sans-serif": "DejaVu Sans", "font.family": "sans-serif"}
)


csv_filename = "Magnify_test_Output_Batch_3d.csv" 
csv_path = Path(csv_filename)

if not csv_path.exists():
    raise FileNotFoundError(f"Soubor '{csv_filename}' nebyl nalezen v aktuální složce!")


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
    raise ValueError(f"Po filtraci (MIN_SIZE={MIN_SIZE}, MAX_SIZE={MAX_SIZE}) nezůstala žádná data kykreslení!")

counts_per_cell = (
    df_filtered.groupby("filename").size().reset_index(name="condensate_count")
)

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

#Graph C: Size vs. Integrated Density
sns.scatterplot(
    data=df_filtered,
    x="size_plot",
    y="integrated_density",
    ax=axes[1, 0],
    color="#2e7d32",
    s=40,
    alpha=0.5,
)
sns.regplot(
    data=df_filtered,
    x="size_plot",
    y="integrated_density",
    ax=axes[1, 0],
    scatter=False,
    color="#1b5e20",
)
axes[1, 0].set_title(title_c, fontsize=13, fontweight="bold")
axes[1, 0].set_xlabel(xlabel_text)
axes[1, 0].set_ylabel("Integrated Density (a.u.)")

#Graph D: Concentration of Signal (Max Intensity)
sns.boxplot(
    y=df_filtered["max_intensity"],
    ax=axes[1, 1],
    color="#f57c00",
    width=0.3,
    boxprops=dict(alpha=0.7),
)
sns.stripplot(
    y=df_filtered["max_intensity"],
    ax=axes[1, 1],
    color="black",
    alpha=0.4,
    jitter=0.2,
    size=4,
)
axes[1, 1].set_title(
    "D) Signal Concentration (Max Intensity)", fontsize=13, fontweight="bold"
)
axes[1, 1].set_ylabel("Max Intensity (a.u.)")

plt.tight_layout()


output_plot = f"{folder_name}_Analysis_Plots.png"
plt.savefig(output_plot, dpi=300)

print(f"Graphs were successfully created and saved as '{output_plot}'")
plt.show()
