# Reports & Visualizations

After processing your images, ExQt automatically generates a four-panel summary figure (`partitioning_plots.py`) to help you explore your condensate data at a glance.

---

## The Four Summary Plots

Only QC-passed (primary) objects are plotted, and K_part only where `K_valid`. Reference lines are working values from `reference_values.py` and are labelled as such.

| Position | Plot | What it shows |
| :--- | :--- | :--- |
| **Top-left** | FA vs. $K_{\text{part}}$ | Whether shape and enrichment co-vary. Shading marks the FA class boundary and $K < 1.5$ (both working values). |
| **Top-right** | Volume vs. $K_{\text{part}}$ | Whether enrichment depends on object size. |
| **Bottom-left** | Volume vs. mean intensity | Raw brightness across sizes, before any ratio. |
| **Bottom-right** | Volume vs. FA | Whether shape depends on size; line at the PSF-like elongation working value. |

## How to Read Each Plot

These are descriptive plots; none of them is by itself a test of a physical model.

### 1. FA vs. K_part
Shows whether more enriched objects tend to be rounder or more elongated. Both axes carry their own artefacts: FA depends on PSF elongation and segmentation (see the FA page), K_part on the nucleoplasm estimate and the detector offset.

### 2. Volume vs. K_part
A flat band means K_part does not depend on size in this data set. Small objects are affected by PSF blurring (their mean intensity mixes with the surroundings), which lowers K_part; a rise of K with size can therefore be optical.

### 3. Volume vs. mean intensity
Raw brightness without normalisation. Useful to spot size-dependent dimming of small objects or saturation of large ones.

### 4. Volume vs. FA
Small objects consist of few voxels and are dominated by the PSF shape, so their FA is less reliable; compare with the PSF-like elongation line and treat trends at small volumes with caution.

## Implementation

- `partitioning_plots.export_partitioning_analysis` (PNG, 2-page PDF, `*_partitioning_summary.json` with an interpretation warning about PSF and Core-Enriched).
- The input must be the QC-annotated frame from `postprocessing._prepare_reporting_frames`; a raw CSV without `primary_qc_valid` is rejected.
