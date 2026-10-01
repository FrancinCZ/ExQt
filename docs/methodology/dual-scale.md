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

## Comparing IF and ExM

ExQt reports sizes on the same biological scale for both methods, so results can be placed side by side. The number and size of detected objects also depend on optical resolution, sample preparation and segmentation; ExQt does not separate these effects, so a difference between the methods is a measured result, not an explanation. Besides the per-object distributions, the total object volume per nucleus (sum over all objects with the same `cell_id`) can be compared.

---

## Implementation

- `Batch.process_condensates` writes `volume_px`, `volume_gel_um3`, `volume_bio_um3`, `equivalent_diameter_um` and the calibration actually applied (`applied_pixel_size_nm`, `applied_z_step_nm`); the source of that calibration is in `*_metadata.json → applied_calibration_by_file`.
- The expansion factor is entered with up to 3 decimals; a rounded ExF changes volumes by ExF³.

## Known Limitations

- **One expansion factor for the whole image.** All lengths are divided by a single, uniform expansion factor. If the nucleus expands differently from the gel, biological sizes are off by that ratio and volumes by its cube. Measure the factor as close as possible to the analysed structures and report its uncertainty. A uniform factor does not change K_part or FA.
- **Axial scale.** The Z-step entered in Settings is used as the true distance between slices. When the refractive index of the sample differs from that of the objective's immersion medium (for example a water-swollen gel imaged with an oil objective), the true axial spacing differs from the nominal Z-step, and ExQt does not correct for this. Objects then appear stretched or compressed along Z, which changes FA and volumes; samples imaged under different conditions (for example IF in mounting medium and ExM in a gel) can differ in this respect. If the axial correction factor for your set-up is known, enter the corrected Z-step in Settings.
