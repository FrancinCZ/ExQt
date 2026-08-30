# 3D Thermodynamic Partitioning ($K_{\text{part}}$) & Biophysical Phase Profiling

This module in **ExQt** provides automated, single-cell quantitative extraction of the **thermodynamic partitioning coefficient ($K_{\text{part}}$)** and generates **2D Biophysical Phase Diagrams** directly from 3D Expansion Microscopy (ExM) confocal stacks.

---

## 1. What is the Partitioning Coefficient ($K_{\text{part}}$)?

The partitioning coefficient quantifies the **thermodynamic propensity of a molecule to enrich within the dense condensed phase relative to the surrounding dilute nucleoplasmic phase**:

$$K_{\text{part}} = \frac{\max(I_{\text{obj}} - I_{\text{offset}}, 0)}{\max(I_{\text{nuc}} - I_{\text{offset}}, 10^{-6})}$$

- **$I_{\text{obj}}$:** Mean or median 3D fluorescence intensity of the segmented condensate.
- **$I_{\text{nuc}}$:** Mean 3D fluorescence intensity of the surrounding background (nucleoplasm or cytosol) within the cell ROI, strictly excluding segmented objects.
- **$I_{\text{offset}}$:** Camera dark noise floor (0.5th percentile of the acquired stack).

> [!NOTE]
> For **$10\times$ Expansion Microscopy ($\text{ExF} = 10$)**, hydrogel volume expands $1000$-fold ($10^3$). Non-crosslinked soluble background is diluted and washed out ($I_{\text{nuc}} \to I_{\text{dark}}$), while densely crosslinked molecules in condensates retain high fluorescence, amplifying the physical signal-to-background contrast.

---

## 2. 2D Biophysical Phase Diagram ($FA$ vs. $K_{\text{part}}$)

By combining **3D Fractional Anisotropy ($FA$)** with **$K_{\text{part}}$**, ExQt categorizes every segmented assembly into distinct biophysical regimes:

```text
       K_part ▲
              │  [QUAD 1: Spherical LLPS]   │  [QUAD 2: Chromatin-Wetted LLPS]
              │  • Unconstrained liquid     │  • Mechanically constrained
              │  • Isotropic surface tension│  • Capillary wetting along DNA
       15-20x ┼── (e.g. SON Speckles) ──────┼── (e.g. POL II Transcription Hubs) ───
              │                             │
              │  [QUAD 3: Weak/Diffuse]     │  [QUAD 4: Hollow/Membrane-Bound]
              │  • Non-condensed background │  • Hollow lumen (GM130 / Golgi)
              │  • No phase separation      │  • Negative layer profile (ΔI < 0)
         1.0x ┼─────────────────────────────┴───────────────────────────────────────
              └─────────────────────────────────────────────────────────────► FA
              0.0 (Perfect Sphere)        0.65               1.0 (Elongated / Sheet)
```


## 3. Output Plots & Deliverables

When partitioning analysis is enabled, ExQt automatically generates:
1. `*_3d_partitioning_analysis.png` — Two-panel figure containing:
   - **Left:** *Size vs. $K_{\text{part}}$* stability plot (proving concentration independence across droplet volumes).
   - **Right:** *2D Biophysical Phase Diagram* with color-coded quadrants.
2. `*_3d_Radial_FA_Profiling_Plots.png` — 4-panel *Core-Middle-Shell* radial profiling and QC Funnel.

---


