import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Nastavení fontu kompatibilního s Linuxem
sns.set_theme(style="ticks", palette="muted")
plt.rcParams.update(
    {"font.sans-serif": "DejaVu Sans", "font.family": "sans-serif"}
)


df = pd.read_csv("Final_Output_Batch.csv")

if "volume_bio_nm3" in df.columns:
    df["size_plot"] = df["volume_bio_nm3"] / 1e9  
    xlabel_text = "Biological Volume (µm³)"
    title_a = "A) Condensate Volume Distribution"
    title_c = "C) Volume vs. Total Intensity"
    unit_label = "µm³"
elif "area_bio_nm2" in df.columns:
    df["size_plot"] = df["area_bio_nm2"] / 1e6 
    xlabel_text = "Biological Area (µm²)"
    title_a = "A) Condensate Area Distribution"
    title_c = "C) Area vs. Total Intensity"
    unit_label = "µm²"
else:
    raise KeyError("V CSV chybí sloupec 'volume_bio_nm3' i 'area_bio_nm2'!")


counts_per_cell = (
    df.groupby("filename").size().reset_index(name="condensate_count")
)


fig, axes = plt.subplots(2, 2, figsize=(15, 11))


sns.histplot(
    df["size_plot"],
    kde=True,
    ax=axes[0, 0],
    color="#2b5c8f",
    bins=30 if len(df) > 30 else 10,
)
median_size = df["size_plot"].median()
axes[0, 0].axvline(
    median_size,
    color="#d9534f",
    linestyle="--",
    linewidth=2,
    label=f"Median: {median_size:.3f} {unit_label}",
)
axes[0, 0].set_title(title_a, fontsize=13, fontweight="bold")
axes[0, 0].set_xlabel(xlabel_text)
axes[0, 0].set_ylabel("Count")
axes[0, 0].legend()

#Grph B
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

#Graph C
sns.scatterplot(
    data=df,
    x="size_plot",
    y="integrated_density",
    ax=axes[1, 0],
    color="#2e7d32",
    s=40,
    alpha=0.5,
)
sns.regplot(
    data=df,
    x="size_plot",
    y="integrated_density",
    ax=axes[1, 0],
    scatter=False,
    color="#1b5e20",
)
axes[1, 0].set_title(title_c, fontsize=13, fontweight="bold")
axes[1, 0].set_xlabel(xlabel_text)
axes[1, 0].set_ylabel("Integrated Density (a.u.)")

#Graph D

sns.boxplot(
    y=df["mean_intensity"],
    ax=axes[1, 1],
    color="#f57c00",
    width=0.3,
    boxprops=dict(alpha=0.7),
)
sns.stripplot(
    y=df["mean_intensity"],
    ax=axes[1, 1],
    color="black",
    alpha=0.4,
    jitter=0.2,
    size=4,
)
axes[1, 1].set_title(
    "D) Signal Concentration (Mean Intensity)", fontsize=13, fontweight="bold"
)
axes[1, 1].set_ylabel("Mean Intensity (a.u.)")

plt.tight_layout()
plt.savefig("ExM_Analysis_Plots.png", dpi=300)
print("Graphs were succesfuly created and saved as'ExM_Analysis_Plots.png'")
plt.show()