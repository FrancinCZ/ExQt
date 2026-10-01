# Z-Stack Aligner (XY Drift Correction)

`stack_aligner.py` corrects **lateral (XY) drift between consecutive Z-slices** of one stack, e.g. from stage or thermal drift during acquisition. It estimates one translation per pair of neighbouring slices from a single reference channel and applies the accumulated correction to **all channels** of the stack (and optionally to the mask). It does not correct chromatic aberration between channels, rotation, scaling or non-rigid deformation.

---

## Algorithm

1. **Plane preparation** (`_prepare_plane`): percentile clipping (1–99 %), high-pass filter (σ = 8 px), light smoothing, normalisation and a 2D Hanning taper against FFT edge effects. Planes without texture are marked `insufficient_texture`.
2. **Registration:** `skimage.registration.phase_cross_correlation` between neighbouring planes, sub-pixel with `upsample_factor = 10` (1/10 px), in both directions; the two estimates must agree (`max_bidirectional_disagreement_px`).
3. **Step checks:** step size ≤ `max_step_px` (default 8 px), correlation after the shift (computed only over overlapping pixels), no loss of correlation. Each step is `PASS`, `REVIEW` or `FAIL`; a `FAIL` step gets zero shift.
4. **Overall status:** texture-less steps at the ends of the stack are ignored; `FAIL` steps in the middle up to `max_fail_fraction` (20 %) give `REVIEW`, more give `FAIL`. The overall status is written to `<stem>_drift.csv` (`overall_status`).
5. **Output:** with `expand_canvas` (default) the canvas is enlarged by the accumulated drift and the new border is filled with zeros. These padding voxels are marked 0 in `<stem>_alignment_valid.tif`, and the analysis excludes them from intensity statistics and uses them for edge detection. For `FAIL` no image is written.

The analysis follows the overall status: `PASS` normal, `REVIEW` analysed but not primary, `FAIL` excluded.

---

## Implementation

- `stack_aligner.estimate_xy_drift`, `apply_xy_shifts`, `align_tiff_stack`, `align_tiff_folder`; the reference channel can be suggested by `ChannelDetectionWorker` in `App.py`.

---

## Methodological Background & Inspiration

The multi-channel 3D stack alignment strategy in `stack_aligner.py` is inspired by and adapted from the MATLAB **[3D-Aligner](https://github.com/suzukilabmcardle/3D-Aligner)** developed by the **Suzuki Lab at the McArdle Laboratory for Cancer Research** ([suzukilabmcardle/3D-Aligner](https://github.com/suzukilabmcardle/3D-Aligner)).

In ExQt, the core concepts were re-engineered and implemented natively in Python using `skimage.registration.phase_cross_correlation` and `scipy.ndimage`, augmented with:
1. **Reference channel suggestion**: **Auto-detect** ranks channels by the correlation between neighbouring slices (`ChannelDetectionWorker`).
2. **Spectral Leakage Suppression**: 2D Hanning window preconditioning to reduce Fourier border discontinuities.
3. **Canvas Expansion (`expand_canvas=True`)**: The canvas grows by the accumulated drift so no signal is cropped; the zero padding is recorded in `<stem>_alignment_valid.tif` and excluded from statistics.
4. **Synchronized Mask Alignment**: Rigid transformation validation and synchronized translation applied automatically to corresponding binary/labeled `*_Mask.tif` stacks.
