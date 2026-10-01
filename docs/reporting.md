# Reports & Visualizations

When reports are enabled, ExQt saves a four-panel summary figure and a PDF report for the primary objects.

---

## The Four Summary Plots

Only QC-passed (primary) objects are plotted, and K_part only where `K_valid`. Reference lines are working values from `reference_values.py` and are labelled as such.

| Position | Plot | What it shows |
| :--- | :--- | :--- |
| **Top-left** | FA vs. $K_{\text{part}}$ | Shading at the FA class boundary and $K = 1.5$ (working values). |
| **Top-right** | Volume vs. $K_{\text{part}}$ | |
| **Bottom-left** | Volume vs. mean intensity | Intensity before any ratio. |
| **Bottom-right** | Volume vs. FA | Line at the PSF-like elongation (working value). |

## How to Read Each Plot

The plots are descriptive; none of them tests a physical model. Keep in mind:

- FA depends on optical blur along Z and on segmentation (see the FA page).
- $K_{\text{part}}$ depends on the nucleoplasm estimate and the detector offset (see the K_part page).
- Small objects consist of few voxels and are blurred, which affects both FA and $K_{\text{part}}$; trends at small volumes can therefore be optical.

## Implementation

- `partitioning_plots.export_partitioning_analysis` (PNG, 2-page PDF, `*_partitioning_summary.json` with an interpretation warning about PSF and Core-Enriched).
- The input must be the QC-annotated frame from `postprocessing._prepare_reporting_frames`; a raw CSV without `primary_qc_valid` is rejected.
