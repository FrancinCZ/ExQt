# Reports & Visualizations

When reports are enabled, ExQt saves a standard six-panel figure (`*_Analysis_Plots.png`), a four-panel summary figure and a PDF report for the primary objects.

---

## The Standard Plots

Every run is drawn with the same six panels and the same fixed axes, so two runs can be compared by eye. Only QC-passed (primary) objects are plotted, and K_part only where `K_valid`. Intensity is shown as $K_{\text{part}}$ (object / nucleoplasm of the same nucleus, baseline 1), not in raw ADU, which is not comparable between images.

| Panel | Plot | Axes |
| :--- | :--- | :--- |
| **A** | Size distribution | Log size axis fixed to the GUI size range. |
| **B** | $K_{\text{part}}$ per nucleus | Log axis fixed to `plot_k_min`–`plot_k_max` (`defaults.py`, working values 0.5–50); line at K = 1. |
| **C** | Size vs. $K_{\text{part}}$ | Same fixed axes as A and B. |
| **D** | FA per nucleus | 0–1; needs Radial FA Profiling. |
| **E** | Size vs. FA | Fixed size axis, FA 0–1. |
| **F** | QC-passed objects per nucleus | With the counts of segmented, size-eligible and QC-passed objects. |

Size axes are logarithmic and labelled in plain decimals (0.01, 0.1, 1 µm³); a range under two decades also gets 0.02 and 0.05 ticks. Points outside a fixed axis are not drawn; the panel states how many there are. The title gives the number of objects and nuclei.

## Conditions Side by Side

**Merge existing runs…** also saves `Merged_Stats_Conditions.png` and a `Per_Nucleus` sheet. A condition is the top-level folder under the chosen merge folder, so replicates in subfolders are grouped:

```
Analysis/
  Control/rep1/…_Output_Batch_3d.csv
  Control/rep2/…
  Treated/rep1/…
```

The figure has the six panels of the standard plots (A size, B $K_{\text{part}}$, C size vs. $K_{\text{part}}$, D FA, E size vs. FA, F objects per nucleus) on the same fixed axes, with one box per condition in A, B, D and F, so several conditions (e.g. Control, Test, Test2) sit side by side in one graph. C and E are scatter plots coloured by condition. In F each dot is one nucleus and the bar is the median. Compare conditions by nuclei: objects of one nucleus are not independent measurements. No statistical test is drawn.

**Tools → Compare runs…** draws the same six panels for groups that cannot be merged, e.g. ExM (expansion factor 8–10) and IF (expansion factor 1). Each chosen folder is one group; runs inside a group must share the QC settings, settings that differ between groups are listed under the figure. Single nuclei can be picked, and each picked nucleus can get its own column (coloured by its group). Unlike the standard plots and Merge, the $K_{\text{part}}$ axis in Compare runs is fitted to the plotted values: from the next round value below the lowest K to the next above the highest (e.g. highest K 1.5 → axis to 2), always with K = 1 inside. All groups in one figure share that axis, but two Compare runs figures can have different K axes. See [GUI overview](user-guide/gui-overview.md).

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
