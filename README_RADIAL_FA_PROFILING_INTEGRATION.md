# ExQt Radial FA Profiling Integration

> Compatibility note: implementation modules, settings, CSV fields, and metadata
> still use the historical internal identifier `mode_a` so older ExQt results remain readable.

## Purpose

This package adds an **optional 3D core-shell Fractional Anisotropy (FA)** analysis mode to ExQt. It measures geometric structural-response metrics of a segmented object. It does **not** directly measure stiffness, liquidity, viscosity, or prove that a segmented region is one biological condensate.

## Files to install

Copy these five files into the same existing ExQt source directory:

```text
app.py
Batch.py
rezim_a_core_shell.py
rezim_a_anisotropy.py
rezim_a_metrics.py
```

Before replacing files, make a backup of the original `app.py` and `Batch.py`. The provided `postprocessing.py` does not need a change for this first integration.

`test_batch_mode_a.py` is optional. It is a developer smoke test and is not required for routine analysis.

For a reproducible numerical/plot test environment with Python 3.12, install `requirements-test.txt` and run:

```text
python test_batch_mode_a.py
```

## GUI workflow

1. Start ExQt normally.
2. Select **Settings → Advanced...**.
3. Enter the correct **Pixel Size XY** and **Z-step** for the selected dataset.
4. Enable **Radial FA Profiling (3D)**.
5. Set **Radial FA Profiling minimum voxels per layer**. The initial pilot value is `20`.
6. Keep **Exclude split Z-slice silhouettes from primary comparison** enabled for conservative primary comparisons.
7. Run ExQt in **3d** mode.

Radial FA Profiling is automatically skipped in `2d` and `single_slice` modes. When it is disabled, the existing ExQt object metrics and CSV structure remain unchanged.

## Calibration rule

The Pixel Size XY and Z-step values saved in ExQt Advanced Settings are the source of truth for Radial FA Profiling. `Batch.py` no longer silently replaces a user-supplied GUI value with embedded TIFF metadata. TIFF metadata is used only when the GUI value is missing.

The same expansion-factor correction used by ExQt biological size metrics is used by Radial FA Profiling:

```text
effective XY = Pixel Size XY / Expansion Factor
effective Z  = Z-step / Expansion Factor
sampling order = (Z, Y, X)
```

## New CSV columns when Radial FA Profiling is enabled

| Column | Meaning |
|---|---|
| `A_object`, `A_shell`, `A_middle`, `A_core` | Fractional Anisotropy of the full object and its radial layers. |
| `Delta_A_middle_shell` | `A_middle - A_shell`; first step of the radial FA profile. |
| `Delta_A_core_middle` | `A_core - A_middle`; second step of the radial FA profile. |
| `Delta_A_core_shell` | `A_core - A_shell`; the full geometric difference, not a mechanical measurement. |
| `A_*_valid` | Whether the FA result passed the minimum voxel QC threshold. |
| `mode_a_core_voxels` | Number of voxels in the core. |
| `mode_a_empty_layers` | Semicolon-separated radial layers that contain no voxels after discretization. Their FA is recorded as `NaN`. |
| `mode_a_layer_complete_coverage` | Whether shell, middle, and core cover every object voxel exactly once. |
| `mode_a_object_touches_edge` | Whether the object touches the image or stack boundary. |
| `mode_a_max_components_per_z_slice` | Largest number of disconnected silhouettes observed in one Z slice. |
| `mode_a_primary_include` | Conservative recommendation for inclusion in a primary biological comparison. |
| `mode_a_qc_reason` | Semicolon-separated QC reasons when primary inclusion is false. |
| `mode_a_sampling_*_nm` | Effective physical sampling used in the actual calculation. |
| `mode_a_layer_scheme` | Versioned definition: `edt_over_max_thirds_v1`. The foreground EDT is divided by its object maximum, then split into thirds. |

The Batch metadata JSON also stores the Radial FA Profiling settings and the interpretation limit.

## QC interpretation

An object is retained in the CSV even when `mode_a_primary_include=False`. This preserves auditability. It should be excluded from the primary comparison if it:

- touches an image or stack edge;
- has more than one disconnected silhouette in any Z slice when the conservative option is enabled; or
- has fewer core voxels than the selected minimum threshold; or
- has an empty shell, middle, or core after discrete radial layering.

The `split_in_z_slice` flag is not proof that two biological condensates were merged. It is a conservative warning that the segmented object may need visual review.

Binary masks may use either 0/1 or 0/255. Masks with several positive integer values are treated as instance-label masks and are rejected if one ID occurs in disconnected 3D components. Probability maps and interpolated, non-integer masks are rejected rather than silently interpreted as object labels.

An empty layer is expected for some very small or thin segmented objects. ExQt writes `NaN` for that layer's FA and a QC reason such as `empty_layer_middle`; it does not stop the entire Batch.

## Validation performed before delivery

`test_batch_mode_a.py` passed with a synthetic 3D TIFF and binary mask. It verified that:

1. Radial FA Profiling disabled adds no `A_*` or radial FA QC columns.
2. Radial FA Profiling enabled writes the expected metrics and QC fields.
3. The edge-object QC flag is recorded.
4. The Batch metric values match the pilot helper on the same local object mask and calibration.
5. The effective sampling follows the ExQt expansion-factor correction.
6. Radial FA Profiling is safely skipped in 2D mode.
7. A thin object with empty radial layers produces `NaN` + QC flags without stopping Batch.

Synthetic data was used only to test software wiring. Biological validation must continue on real masks with visual QC and appropriately designed experimental controls.

When readable OME/ImageJ/TIFF calibration metadata are present and consistent across the batch, ExQt offers to place them into Advanced Settings. It never applies them silently, and the effective XY/Z sampling is shown for explicit confirmation before every run.
