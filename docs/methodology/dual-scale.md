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

## Implementation

- `Batch.process_condensates` writes `volume_px`, `volume_gel_um3`, `volume_bio_um3`, `equivalent_diameter_um` and the calibration actually applied (`applied_pixel_size_nm`, `applied_z_step_nm`); the source of that calibration is in `*_metadata.json → applied_calibration_by_file`.
- The expansion factor is entered with up to 3 decimals; a rounded ExF changes volumes by ExF³.

## Known Limitations

- **Gel ExF vs nuclear ExF.** A single global ExF is applied. Nuclei can expand less than the gel (Pesce et al. 2019, *J. Biophotonics*: gel 5.0×, NPC radius 4.3×, intra-nuclear pore distances 3.8×). Use an ExF measured on nuclei where possible and report its uncertainty; volumes scale with ExF³. K_part and FA do not depend on a uniform ExF.
- **Axial scaling (refractive-index mismatch).** A water-swollen gel (n ≈ 1.33) imaged with an oil objective (n = 1.515) has a true axial step of ≈ 0.85 × the nominal Z-step (Besseling et al. 2015, *J. Microsc.*). ExQt does not yet correct for this: objects appear ~18 % elongated in Z (a sphere gives FA ≈ 0.2) and volumes are overestimated. IF samples in mounting medium have a different factor, so FA and volumes are not directly comparable between IF and ExM without correction.
