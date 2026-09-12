# Dual-Scale Resolution & Volume Normalization

In Expansion Microscopy (ExM), biological specimens are physically magnified within a swellable hydrogel polymer network. Consequently, two distinct spatial frames of reference exist simultaneously:

1. **Hydrogel Space (Physical Measurement Frame):** The physical coordinates recorded by the microscope camera post-expansion.
2. **Biological Space (Pre-Expansion Reference Frame):** The estimated native cellular scale prior to hydrogel anchoring and swelling.

To avoid reporting ambiguity, ExQt natively preserves and exports both scales for every detected object.

---

## Mathematical Formulation

Let $\Delta x, \Delta y$ denote the lateral pixel dimensions (in nm) and $\Delta z$ the axial optical section spacing (in nm) recorded by the microscope. For an isotropic linear expansion factor $\text{ExF}$, the effective biological sampling increments are:


$$\Delta x_{\text{bio}} = \frac{\Delta x}{\text{ExF}}, \quad \Delta y_{\text{bio}} = \frac{\Delta y}{\text{ExF}}, \quad \Delta z_{\text{bio}} = \frac{\Delta z}{\text{ExF}}$$

For a segmented 3D condensate comprising $N_{\text{voxels}}$:

$$\text{Vol}_{\text{gel}} = N_{\text{voxels}} \times \frac{\Delta x \cdot \Delta y \cdot \Delta z}{10^9} \quad [\mu\text{m}^3]$$

$$\text{Vol}_{\text{bio}} = N_{\text{voxels}} \times \frac{\Delta x_{\text{bio}} \cdot \Delta y_{\text{bio}} \cdot \Delta z_{\text{bio}}}{10^9} = \frac{\text{Vol}_{\text{gel}}}{\text{ExF}^3} \quad [\mu\text{m}^3]$$

The equivalent spherical diameter ($D_{\text{eq}}$) in biological space is:

$$D_{\text{eq}} = \left( \frac{6 \cdot \text{Vol}_{\text{bio}}}{\pi} \right)^{1/3} \quad [\mu\text{m}]$$

---

## Comparison Between Imaging Modalities

A common empirical observation when comparing classical IF to ExM is a shift in the median detected object volume: objects that appear as single large entities at diffraction-limited resolution can resolve into several smaller, discrete objects after expansion.

This shift is consistent with two non-exclusive interpretations:

1. **Sub-domain resolution** — ExM resolves constituent sub-compartments that were previously merged within the PSF envelope.
2. **Structural artifacts** — differences in fixation, antibody penetration, or gel homogeneity across conditions may independently affect apparent object sizes.

As a consistency check, total nuclear condensate volume (sum across all detected objects per nucleus) can be compared between modalities. Agreement between IF and ExM totals is consistent with — but does not prove — the sub-domain hypothesis; it primarily indicates that the segmentation and calibration are internally consistent across protocols.



---

## Production Code: Dual-Scale Implementation (`Batch.py`)

Below is the production implementation from [`Batch.py:621–635`](file:///C:/Users/franc/Desktop/ExQt_Rezim_A_Final/Batch.py#L621-L635):

```python
#EXCERPT FROM Batch.py: process_batch_pair()

#Convert pixel counts into calibrated biological units while keeping raw counts for auditability
row["applied_pixel_size_nm"] = float(pixel_size_nm)
row["applied_z_step_nm"] = float(z_step_nm) if is_3d else float("nan")

if is_3d:
    voxel_volume_bio_um3 = ((eff_pixel_size_nm ** 2) * eff_z_step_nm) / 1e9
    voxel_volume_gel_um3 = ((float(pixel_size_nm) ** 2) * float(z_step_nm)) / 1e9
    
    row["volume_px"] = region.area
    row["volume_bio_um3"] = float(region.area * voxel_volume_bio_um3)
    row["volume_gel_um3"] = float(region.area * voxel_volume_gel_um3)
    row["shape_metric_bio"] = row["volume_bio_um3"]
    row["equivalent_diameter_um"] = float((6.0 * row["volume_bio_um3"] / np.pi) ** (1.0 / 3.0))
else:
    pixel_area_bio_um2 = (eff_pixel_size_nm ** 2) / 1e6
    pixel_area_gel_um2 = (float(pixel_size_nm) ** 2) / 1e6
    
    row["area_px"] = region.area
    row["area_bio_um2"] = float(region.area * pixel_area_bio_um2)
    row["area_gel_um2"] = float(region.area * pixel_area_gel_um2)
    row["shape_metric_bio"] = row["area_bio_um2"]
```