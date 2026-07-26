import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


sns.set_theme(style="ticks", palette="muted")
plt.rcParams.update({"font.sans-serif": "Arial", "font.family": "sans-serif"})


df = pd.read_csv("Final_Output_Batch.csv")


if "volume_bio_nm3" in df.columns:
    df["volume_bio_um3"] = df["volume_bio_nm3"] / 1e9
else:

    df["volume_bio_um3"] = df["area_bio_nm2"] / 1e6


counts_per_cell = df.groupby("filename").size().reset_index(name="condensate_count")


fig, axes = plt.subplots(2, 2, figsize=(15, 11))

sns.histplot(
    df["volume_bio_um3"],
    kde=True,
    ax=axes[0, 0],
    color="#2b5c8f",
    bins=30 if len(df) > 30 else 10,
)
median_vol = df["volume_bio_um3"].median()
axes[0, 0].axvline(
    median_vol,
    color="#d9534f",
    linestyle="--",
    linewidth=2,
    label=f"Median: {median_vol:.3f} µm³",
)
axes[0, 0].set_title("A) Condensate Volume Distribution", fontsize=13, fontweight="bold")
axes[0, 0].set_xlabel("Biological Volume (µm³)")
axes[0, 0].set_ylabel("Count")
axes[0, 0].legend()


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
axes[0, 1].set_title("B) Number of Condensates per Cell (FOV)", fontsize=13, fontweight="bold")
axes[0, 1].set_ylabel("Number of Condensates")


sns.scatterplot(
    data=df,
    x="volume_bio_um3",
    y="integrated_density",
    ax=axes[1, 0],
    color="#2e7d32",
    s=40,
    alpha=0.5,
)
sns.regplot(
    data=df,
    x="volume_bio_um3",
    y="integrated_density",
    ax=axes[1, 0],
    scatter=False,
    color="#1b5e20",
)
axes[1, 0].set_title("C) Volume vs. Total Intensity", fontsize=13, fontweight="bold")
axes[1, 0].set_xlabel("Biological Volume (µm³)")
axes[1, 0].set_ylabel("Integrated Density (a.u.)")


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
axes[1, 1].set_title("D) Signal Concentration (Mean Intensity)", fontsize=13, fontweight="bold")
axes[1, 1].set_ylabel("Mean Intensity (a.u.)")

plt.tight_layout()
plt.savefig("ExM_3D_Analysis_Plots.png", dpi=300)
print("Updated plots saved successfully as 'ExM_3D_Analysis_Plots.png'")
plt.show()