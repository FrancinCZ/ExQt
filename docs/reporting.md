# Reports & Visualizations

After processing your images, ExQt automatically generates a four-panel summary figure (`partitioning_plots.py`) to help you explore your condensate data at a glance.

---

## The Four Summary Plots

The generated dashboard contains four key plots:

| Position | What it plots | Why it matters in plain terms |
| :--- | :--- | :--- |
| **Top-Left** | **Shape vs. Enrichment**<br/>(FA vs. $K_{\text{part}}$) | Do denser droplets have different shapes? Checks if protein concentration relates to whether the condensate is round or elongated. |
| **Top-Right** | **Size vs. Enrichment**<br/>(Volume vs. $K_{\text{part}}$) | **The classic droplet test:** In a true liquid droplet, concentration should stay roughly constant no matter how big or small the droplet grows. |
| **Bottom-Left** | **Volume vs. Raw Brightness**<br/>(Volume vs. Mean Intensity) | **Sanity check:** Verifies the camera signal across sizes to ensure measurements aren't biased by microscope settings. |
| **Bottom-Right** | **Size vs. Shape**<br/>(Volume vs. FA) | **Chromatin confinement:** Small droplets are usually round (like tiny raindrops). Larger droplets often get squeezed and deformed by surrounding chromatin. |

---

## How to Read Each Plot

### 1. Shape vs. Enrichment (FA vs. K_part)
- **What is FA and $K_{\text{part}}$?** Fractional Anisotropy ($\text{FA}$) measures elongation (0 = sphere, higher = stretched out), while $K_{\text{part}}$ measures protein enrichment inside the droplet relative to the nucleoplasm.
- **What to look for:** If droplets behave like simple liquids, surface tension keeps them relatively round ($\text{FA} < 0.4$). If they form more rigid, elongated scaffolds (like nuclear speckles), you may see higher anisotropy values.

### 2. Size vs. Enrichment (Volume vs. K_part)
- **The core idea:** Think of vinegar and oil. A large droplet of oil has the same concentration ($K_{\text{part}}$) of oil molecules as a small droplet.
- **What to look for:** In pure liquid-like condensates, data points form a flat, horizontal band across different volumes.

### 3. Volume vs. Raw Intensity
- **The core idea:** This looks at raw camera brightness before calculating ratios.
- **What to look for:** Ensures your signal is consistent across different sizes and not artificially dimming in smaller objects.

### 4. Size vs. Shape (Volume vs. FA)
- **The core idea:** Inside the nucleus, condensates are embedded in a mesh of chromatin fibers.
- **What to look for:** 
  - **Small condensates** have enough surface tension to easily push chromatin aside and stay round (low FA).
  - **Larger condensates** cannot push chromatin away as easily, so they grow into gaps and become irregular or elongated (higher FA).


---

## Production Code: Four-Panel Figure Exporter (`partitioning_plots.py`)

Below is the production implementation from [`partitioning_plots.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/partitioning_plots.py) that generates the four-panel figure:

```python
#EXCERPT FROM partitioning_plots.py: export_partitioning_report()

def export_partitioning_report(
    df: pd.DataFrame,
    output_png: Path | str,
    title: str = "Condensate Quantitative Analysis",
    dpi: int = 300,
) -> Path:
    """Generates a unified 2x2 multi-panel publication dashboard."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), dpi=dpi)
    
    #1. Morphology vs. Partitioning (Top-Left)
    plot_morphology_vs_partitioning(df, ax=axes[0, 0])
    
    #2. Size vs. Partitioning (Top-Right — Thermodynamic LLPS Test)
    plot_size_vs_partitioning(df, ax=axes[0, 1])
    
    #3. Volume vs. Mean Fluorescence Intensity (Bottom-Left)
    plot_volume_vs_intensity(df, ax=axes[1, 0])
    
    #4. Volume vs. 3D Shape Anisotropy (Bottom-Right — Chromatin Confinement)
    plot_volume_vs_anisotropy(df, ax=axes[1, 1])
    
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    out_path = Path(output_png)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path
```