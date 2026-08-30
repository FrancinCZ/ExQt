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

### Four Biophysical Regimes:
1. **Spherical LLPS Droplets ($FA < 0.65, K_{\text{part}} \ge 1.5$):** Unconstrained liquid droplets dominated by isotropic surface tension (e.g., **SON Nuclear Speckles**, $FA \approx 0.533, K_{\text{part}} \approx 18.57\times$).
2. **Chromatin-Wetted Condensates ($FA \ge 0.65, K_{\text{part}} \ge 1.5$):** Condensed assemblies nucleating on linear DNA polymers, exhibiting capillary elongation along the chromatin scaffold (e.g., **POL II Hubs**, $FA \approx 0.810, K_{\text{part}} \approx 18.93\times$).
3. **Weak / Diffuse Assemblies ($K_{\text{part}} < 1.5$):** Structures failing to maintain dense phase concentration above background.
4. **Hollow / Membrane-Bound Organelles:** Hollow 2D membranes enclosing empty non-fluorescent lumens (e.g., **Golgi apparatus**, $K_{\text{part}} = 0.89\times, 31.1\%$ hollow layers $\Delta I < 0$).

---

## 3. Experimental Benchmarks

| Metric (3D Analysis) | 🔵 POL II (Transcriptional Hubs) | 🔴 GOLGI / GM130 (Negative Control) | 🟢 SON (Nuclear Speckles) |
| :--- | :---: | :---: | :---: |
| **Biological State** | Chromatin-tethered condensate | Hollow membrane labyrinth | Unconstrained spherical droplet |
| **Median Intensity ($I_{\text{obj}}$)** | **$4\,724.72\text{ ADU}$** | **$679.46\text{ ADU}$** | **$5\,120.30\text{ ADU}$** |
| **Background ($I_{\text{nuc}}$)** | **$249.58\text{ ADU}$** | **$761.23\text{ ADU}$** | **$275.40\text{ ADU}$** |
| **Partitioning ($K_{\text{part}}$)** | **$\mathbf{18.93\times}$** | **$\mathbf{0.89\times}$** | **$\mathbf{18.57\times}$** |
| **Fractional Anisotropy ($FA$)** | **$0.810$** *(prolate ellipsoid)* | **$0.698\text{ to }0.99$** *(planar sheets)* | **$0.533$** *(globular spheroid)* |
| **Hollow Structures ($\Delta I < 0$)** | **$0\,\%$** | **$31.1\,\%$ in large cisternae** | **$0\,\%$** |

---

## 4. Output Plots & Deliverables

When partitioning analysis is enabled, ExQt automatically generates:
1. `*_3d_partitioning_analysis.png` — Two-panel figure containing:
   - **Left:** *Size vs. $K_{\text{part}}$* stability plot (proving concentration independence across droplet volumes).
   - **Right:** *2D Biophysical Phase Diagram* with color-coded quadrants.
2. `*_3d_Radial_FA_Profiling_Plots.png` — 4-panel *Core-Middle-Shell* radial profiling and QC Funnel.
3. Full validation reports (`exm_partitioning_report_EN.docx` / `.html` / `.md`).

---

## 5. References

1. **Cho, Spille, Cisse et al. (*Science* 2018):** Quantification of 200–400 molecules in $>300\text{ nm}$ transcription-dependent condensates ([PMC6543815](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6543815/)).
2. **Feric, Pappu, Brangwynne et al. (*Cell* 2016):** Thermodynamic hierarchy and surface tension driving core-shell multiphase condensate architecture ([DOI: 10.1016/j.cell.2016.04.047](https://doi.org/10.1016/j.cell.2016.04.047)).
3. **Wang, Zhang et al. (*Biophysics Reports* 2018):** Standardized protocol for quantitative protein phase partitioning and wetting assays ([DOI: 10.1007/s41048-018-0078-7](https://doi.org/10.1007/s41048-018-0078-7)).
