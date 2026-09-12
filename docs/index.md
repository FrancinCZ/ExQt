# ExQt: Expansion Microscopy Quantification Tool

**ExQt** (*Expansion Microscopy Quantification Tool*) is an open-source, automated image-analysis platform for the 3D segmentation, biophysical quantification, and morphological profiling of biomolecular condensates.

It is designed to work across both **classical confocal immunofluorescence (IF)** and **Expansion Microscopy (ExM)**, providing a unified, auditable pipeline to investigate liquid–liquid phase separation, internal compartmentalization, and morphological organization in the mammalian cell nucleus.

---

## Why ExQt?

In diffraction-limited confocal microscopy, structures closer than the Abbe resolution limit can overlap and appear as a single object. Expansion Microscopy addresses this by physically magnifying the specimen in a swellable polymer hydrogel, allowing standard confocal hardware to resolve finer structural detail.

Analyzing ExM data introduces specific computational requirements:

- **True 3D processing** — maintaining connectivity across anisotropic voxel grids ($\Delta z > \Delta x, \Delta y$).
- **Dual-scale reporting** — preserving both the physical hydrogel scale ($V_{\text{gel}}$) and the estimated biological pre-expansion scale ($V_{\text{bio}} = V_{\text{gel}}/\text{ExF}^3$).
- **Internal profiling** — quantifying concentration gradients across concentric layers (Core–Middle–Shell), rather than treating condensates as uniform blobs.
- **Shape analysis** — separating external elongation from genuine internal organization via geometric reference models.

---

## Pipeline Orchestration (`Batch.py`)

The analysis kernel is `process_batch_pair()` in [`Batch.py`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/Batch.py). The excerpt below shows how voxel counts are converted into calibrated biophysical measurements:

```python
#Batch.py: process_batch_pair()

#Effective biological pixel scale (after expansion correction)
eff_pixel_size_nm = pixel_size_nm / expansion_factor
eff_z_step_nm     = z_step_nm / expansion_factor

#Camera dark offset: lowest 0.5th percentile inside the nuclear ROI
roi_pixels    = img_intensity[roi_mask] if np.any(roi_mask) else img_intensity.ravel()
camera_offset = float(np.percentile(roi_pixels, 0.5)) if roi_pixels.size > 0 else 0.0

#Nuclear background (nucleoplasm: ROI interior minus condensate voxels)
nuc_pixels       = img_intensity[roi_mask & (~(labeled_mask > 0))]
nucleoplasm_mean = float(np.mean(nuc_pixels)) if nuc_pixels.size > 0 else float("nan")

#Per-object metric loop
for region in measure.regionprops(labeled_mask, intensity_image=img_intensity):
    voxel_volume_bio_um3 = ((eff_pixel_size_nm ** 2) * eff_z_step_nm) / 1e9
    voxel_volume_gel_um3 = ((float(pixel_size_nm) ** 2) * float(z_step_nm)) / 1e9

    row["volume_px"]              = region.area
    row["volume_bio_um3"]         = float(region.area * voxel_volume_bio_um3)
    row["volume_gel_um3"]         = float(region.area * voxel_volume_gel_um3)
    row["equivalent_diameter_um"] = float((6.0 * row["volume_bio_um3"] / np.pi) ** (1.0 / 3.0))

    #Partition coefficient (NaN if denominator is near-zero)
    denom = nucleoplasm_mean - camera_offset
    if np.isnan(nucleoplasm_mean) or denom <= 1.0:
        row["partition_coefficient"] = float("nan")
        row["K_valid"] = False
    else:
        k_val = (region.mean_intensity - camera_offset) / denom
        row["partition_coefficient"] = float(k_val) if k_val >= 1.0 else float("nan")
        row["K_valid"] = k_val >= 1.0

    #Mode A: radial profiling and geometric null model
    if mode_a_enabled:
        radial_metrics = compute_core_shell_metrics(
            region.image.astype(bool),
            region.intensity_image,
            voxel_spacing=(eff_z_step_nm, eff_pixel_size_nm, eff_pixel_size_nm)
        )
        row.update(radial_metrics)
```

---

## Installation

### Prerequisites
- Python 3.10–3.14
- Windows, Linux, or macOS

### Setup
```bash
git clone https://github.com/FrancinCZ/ExQt.git
cd ExQt
pip install -r requirements.txt
python app.py
```

---
