# Invariants of Batch.process_condensates and the geometry helpers on synthetic data with a
# known answer, plus regression tests for defects found in review.
import os
import warnings

import numpy as np
import pandas as pd
import pytest
import tifffile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from Batch import (  # noqa: E402
    _physical_value_to_nm,
    _prepare_labeled_mask,
    _touches_roi_edge,
    get_metadata_from_tif,
    nucleus_roi_labels,
    process_condensates,
)
from skimage.measure import label, regionprops  # noqa: E402
from rezim_a_anisotropy import shape_anisotropy  # noqa: E402

XY_NM, Z_NM = 100.0, 300.0


def _write(folder, raw, mask, name="img", raw_axes="ZYX", mask_axes="ZYX"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        tifffile.imwrite(folder / f"{name}.tif", raw, metadata={"axes": raw_axes})
        tifffile.imwrite(folder / f"{name}_Mask.tif", mask, metadata={"axes": mask_axes})
    return folder / f"{name}.tif", folder / f"{name}_Mask.tif"


def _run(tmp_path, raw, mask, roi=None, *, mode="3d", raw_axes="ZYX", mask_axes="ZYX", **kwargs):
    tif, msk = _write(tmp_path, raw, mask, raw_axes=raw_axes, mask_axes=mask_axes)
    options = dict(
        mode=mode, min_voxels=1, pixel_size_nm=XY_NM, z_step_nm=Z_NM,
        signal_channel=0, dapi_channel=0, roi_output_dir=tmp_path / "out", show_napari=False,
    )
    if roi is not None:
        options.update(auto_roi=False, request_roi_func=lambda shape, is_3d: roi.copy())
    options.update(kwargs)
    return process_condensates(tif, msk, **options)


def _one_object_scene(object_value=700, background=100):
    raw = np.full((8, 40, 40), background, dtype=np.uint16)
    mask = np.zeros(raw.shape, dtype=np.uint8)
    mask[3:5, 15:20, 15:20] = 1
    raw[mask > 0] = object_value
    roi = np.zeros(raw.shape, dtype=np.int64)
    roi[4, 5:35, 5:35] = 1                        # one slice -> extruded over all Z
    return raw, mask, roi


# --- K_part ---------------------------------------------------------------------------------

def test_k_part_with_offset_is_exact(tmp_path):
    raw, mask, roi = _one_object_scene(700, 100)
    df = _run(tmp_path, raw, mask, roi, detector_offset_adu=50.0)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["K_valid"] and row["K_invalid_reason"] == ""
    assert row["nucleoplasm_mean_intensity"] == pytest.approx(100.0)
    assert row["partition_coefficient"] == pytest.approx((700 - 50) / (100 - 50))


def test_offset_at_or_above_nucleoplasm_gives_no_k(tmp_path):
    raw, mask, roi = _one_object_scene(700, 100)
    df = _run(tmp_path, raw, mask, roi, detector_offset_adu=100.0)
    row = df.iloc[0]
    assert not row["K_valid"]
    assert row["K_invalid_reason"] == "negative_denominator"
    assert np.isnan(row["partition_coefficient"])


def test_k_is_independent_of_global_intensity_scale_without_offset(tmp_path):
    raw, mask, roi = _one_object_scene(700, 100)
    k1 = _run(tmp_path / "a", *_mkdirs(tmp_path / "a", raw, mask, roi)).iloc[0]["partition_coefficient"]
    k3 = _run(tmp_path / "b", *_mkdirs(tmp_path / "b", raw * 3, mask, roi)).iloc[0]["partition_coefficient"]
    assert k1 == pytest.approx(7.0) and k3 == pytest.approx(k1)


def _mkdirs(folder, raw, mask, roi):
    folder.mkdir(parents=True, exist_ok=True)
    return raw, mask, roi


# --- channel selection and axis order --------------------------------------------------------

@pytest.mark.parametrize("axes", ["ZCYX", "CZYX"])
def test_signal_channel_is_taken_from_the_c_axis_in_any_order(tmp_path, axes):
    raw, mask, roi = _one_object_scene(700, 100)
    decoy = np.full_like(raw, 5000)               # channel 0: would give K = 1
    stack = np.stack([decoy, raw], axis=axes.index("C"))
    df = _run(tmp_path, stack, mask, roi, raw_axes=axes, signal_channel=1, dapi_channel=0)
    assert df.iloc[0]["partition_coefficient"] == pytest.approx(7.0)


def test_out_of_range_signal_channel_is_an_error(tmp_path):
    raw, mask, roi = _one_object_scene()
    stack = np.stack([raw, raw], axis=1)
    with pytest.raises(ValueError, match="Signal channel"):
        _run(tmp_path, stack, mask, roi, raw_axes="ZCYX", signal_channel=2)


# --- sizes and expansion factor ---------------------------------------------------------------

def test_biological_volume_is_gel_volume_over_exf_cubed(tmp_path):
    raw, mask, roi = _one_object_scene()
    df = _run(tmp_path, raw, mask, roi, expansion_factor=4.0)
    row = df.iloc[0]
    voxels = int(mask.sum())
    assert row["volume_px"] == voxels
    assert row["volume_gel_um3"] == pytest.approx(voxels * XY_NM * XY_NM * Z_NM / 1e9)
    assert row["volume_bio_um3"] == pytest.approx(row["volume_gel_um3"] / 4.0**3)


# --- reproducibility (VALIDATION_PROTOCOL §2) --------------------------------------------------

def test_repeated_runs_give_identical_tables(tmp_path):
    rng = np.random.default_rng(7)
    raw = rng.normal(200, 20, (12, 48, 48)).clip(0).astype(np.uint16)
    zz, yy, xx = np.mgrid[:12, :48, :48]
    blob = ((zz - 6) * 3.0) ** 2 + (yy - 24) ** 2 / 1.5 + (xx - 24) ** 2 <= 49
    mask = blob.astype(np.uint8)
    raw[blob] += 600
    roi = np.zeros(raw.shape, dtype=np.int64)
    roi[6, 4:44, 4:44] = 1
    first = _run(tmp_path / "r1", *_mkdirs(tmp_path / "r1", raw, mask, roi), mode_a_enabled=True)
    second = _run(tmp_path / "r2", *_mkdirs(tmp_path / "r2", raw, mask, roi), mode_a_enabled=True)
    assert len(first) == 1 and "A_core" in first.columns
    pd.testing.assert_frame_equal(first, second)


# --- mask validation -----------------------------------------------------------------------

def test_binary_0_255_mask_is_split_into_connected_objects():
    mask = np.zeros((4, 20, 20), dtype=np.uint8)
    mask[1:3, 2:5, 2:5] = 255
    mask[1:3, 10:14, 10:14] = 255
    assert int(_prepare_labeled_mask(mask).max()) == 2


def test_reused_instance_label_on_disconnected_objects_is_an_error():
    mask = np.zeros((4, 20, 20), dtype=np.uint16)
    mask[1:3, 2:5, 2:5] = 1
    mask[1:3, 10:14, 10:14] = 1
    mask[1:3, 16:19, 16:19] = 2
    with pytest.raises(ValueError, match="disconnected"):
        _prepare_labeled_mask(mask)


def test_probability_mask_is_rejected():
    with pytest.raises(ValueError, match="non-integer"):
        _prepare_labeled_mask(np.full((2, 4, 4), 0.5))


# --- units ------------------------------------------------------------------------------------

@pytest.mark.parametrize("value, unit, expected", [
    (0.1, "µm", 100.0), (0.1, "μm", 100.0), (0.1, "micron", 100.0), (0.1, " UM ", 100.0),
    (58, "nm", 58.0), (0.0001, "mm", 100.0),
])
def test_physical_units_are_converted_to_nm(value, unit, expected):
    assert _physical_value_to_nm(value, unit) == pytest.approx(expected)


@pytest.mark.parametrize("value, unit", [(0.1, "pixel"), (0.1, ""), (0.0, "um"), (-1.0, "um"), (1.0, None)])
def test_unknown_unit_or_non_positive_value_is_not_a_calibration(value, unit):
    assert _physical_value_to_nm(value, unit) is None


# --- shape anisotropy -------------------------------------------------------------------------

def _ellipsoid(semi_axes, angle_deg=0.0, size=61):
    c = size // 2
    zz, yy, xx = np.mgrid[:size, :size, :size] - c
    a = np.deg2rad(angle_deg)
    u = np.cos(a) * xx + np.sin(a) * yy           # rotation in the XY plane
    v = -np.sin(a) * xx + np.cos(a) * yy
    return (zz / semi_axes[0]) ** 2 + (v / semi_axes[1]) ** 2 + (u / semi_axes[2]) ** 2 <= 1.0


@pytest.mark.parametrize("angle", [0.0, 30.0, 45.0, 77.0])
def test_fa_is_rotation_invariant(angle):
    reference = shape_anisotropy(_ellipsoid((8, 8, 20)))["fractional_anisotropy"]
    rotated = shape_anisotropy(_ellipsoid((8, 8, 20), angle))["fractional_anisotropy"]
    assert rotated == pytest.approx(reference, abs=0.01)


def test_globular_boundary_geometry_matches_its_description():
    # reference_values: "FA 0.65 ≈ axis ratio 1.84 of a prolate ellipsoid"
    eig = np.array([1.84**2, 1.0, 1.0])
    fa = np.sqrt(1.5 * np.sum((eig - eig.mean()) ** 2) / np.sum(eig**2))
    assert fa == pytest.approx(0.65, abs=0.005)


# --- regressions ---------------------------------------------------------------------------------

@pytest.mark.parametrize("n_z", [3, 4])
def test_saved_roi_with_few_z_planes_is_a_grayscale_stack(tmp_path, n_z):
    raw = np.full((n_z, 30, 30), 100, dtype=np.uint16)
    mask = np.zeros(raw.shape, dtype=np.uint8)
    mask[1, 10:14, 10:14] = 1
    roi = np.zeros(raw.shape, dtype=np.int64)
    roi[1, 5:25, 5:25] = 1
    _run(tmp_path, raw, mask, roi)
    with tifffile.TiffFile(tmp_path / "out" / "img_ROI.tif") as tif:
        assert len(tif.pages) == n_z
        assert tif.pages[0].photometric == tifffile.PHOTOMETRIC.MINISBLACK


def _touching_2d_rois():
    roi = np.zeros((40, 40), dtype=np.int64)
    roi[2:38, 2:20] = 1
    roi[2:38, 20:38] = 2
    return roi


def test_touching_2d_nucleus_rois_are_rejected_like_in_3d(tmp_path):
    raw = np.full((6, 40, 40), 100, dtype=np.uint16)
    mask = np.zeros(raw.shape, dtype=np.uint8)
    mask[2:4, 15:19, 17:23] = 1
    with pytest.raises(ValueError, match="overlap or touch"):
        _run(tmp_path, raw, mask, _touching_2d_rois(), mode="2d")


def test_separate_2d_shapes_become_nuclei_in_reading_order():
    drawn = np.zeros((40, 40), dtype=np.int64)
    drawn[25:35, 5:15] = 1                        # drawn first, but lower in the image
    drawn[5:15, 25:35] = 2
    cells = nucleus_roi_labels(drawn)
    assert cells[30, 10] == 2 and cells[10, 30] == 1
    assert set(np.unique(cells)) == {0, 1, 2}


def test_object_reaching_into_a_neighbouring_nucleus_roi_touches_its_roi_edge():
    cells = _touching_2d_rois()                   # e.g. an older saved label image
    objects = np.zeros((40, 40), dtype=np.uint8)
    objects[15:19, 17:23] = 1
    region = regionprops(label(objects))[0]
    assert not _touches_roi_edge(region, cells > 0)   # the union ROI does not see the cut
    assert _touches_roi_edge(region, cells, cell_id=1)
    assert _touches_roi_edge(region, cells, cell_id=2)
    inside = np.zeros((40, 40), dtype=np.uint8)
    inside[15:19, 5:10] = 1
    assert not _touches_roi_edge(regionprops(label(inside))[0], cells, cell_id=1)


def test_object_across_two_nuclei_gets_the_majority_nucleus(tmp_path):
    raw = np.full((6, 40, 40), 100, dtype=np.uint16)
    mask = np.zeros(raw.shape, dtype=np.uint8)
    mask[2:4, 10:20, 10:20] = 1                   # diagonal contact of two blocks = one object
    mask[2:4, 20:22, 20:22] = 1
    raw[mask > 0] = 500
    roi = np.zeros((6, 40, 40), dtype=np.int64)
    roi[3, 2:20, 2:20] = 1                        # nuclei touch only at a corner
    roi[3, 20:38, 20:38] = 2
    df = _run(tmp_path, raw, mask, roi)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["cell_id"] == 1 and bool(row["touches_roi_edge"])


def test_concave_object_with_centroid_outside_roi_is_not_silently_dropped(tmp_path):
    raw = np.full((6, 40, 40), 100, dtype=np.uint16)
    roi = np.zeros(raw.shape, dtype=np.int64)
    roi[3, 5:35, 5:35] = 1
    roi[3, 12:35, 15:25] = 0                      # U-shaped nucleus outline
    mask = np.zeros(raw.shape, dtype=np.uint8)
    mask[2:4, 8:30, 10:14] = 1
    mask[2:4, 8:30, 26:30] = 1
    mask[2:4, 8:12, 10:30] = 1                    # U-shaped object around the notch
    raw[mask > 0] = 500
    df = _run(tmp_path, raw, mask, roi)
    assert len(df) == 1


@pytest.mark.parametrize("mode", ["3D", "mip"])
def test_unknown_mode_is_a_clear_error(tmp_path, mode):
    raw, mask, roi = _one_object_scene()
    with pytest.raises(ValueError, match="mode"):
        _run(tmp_path, raw, mask, roi, mode=mode)


def test_non_square_pixels_are_not_reported_as_one_xy_size(tmp_path):
    raw = np.zeros((4, 16, 16), dtype=np.uint16)
    tifffile.imwrite(tmp_path / "ns.ome.tif", raw, ome=True, metadata={
        "axes": "ZYX", "PhysicalSizeX": 0.1, "PhysicalSizeY": 0.2, "PhysicalSizeZ": 0.3,
    })
    meta = get_metadata_from_tif(tmp_path / "ns.ome.tif") or {}
    assert "pixel_size" not in meta
    assert meta["non_square_pixels"] == (pytest.approx(100.0), pytest.approx(200.0))


def test_non_square_pixels_exclude_the_file_even_with_gui_calibration(tmp_path):
    raw = np.full((4, 16, 16), 100, dtype=np.uint16)
    mask = np.zeros(raw.shape, dtype=np.uint8)
    mask[1:3, 5:9, 5:9] = 1
    tifffile.imwrite(tmp_path / "ns.ome.tif", raw, ome=True, metadata={
        "axes": "ZYX", "PhysicalSizeX": 0.1, "PhysicalSizeY": 0.2, "PhysicalSizeZ": 0.3,
    })
    tifffile.imwrite(tmp_path / "ns_Mask.tif", mask)
    with pytest.raises(ValueError, match="non-square pixels"):
        process_condensates(tmp_path / "ns.ome.tif", tmp_path / "ns_Mask.tif", mode="3d",
                            pixel_size_nm=XY_NM, z_step_nm=Z_NM, signal_channel=0, dapi_channel=0,
                            roi_output_dir=tmp_path / "out")


def test_square_pixels_with_rounding_noise_are_accepted(tmp_path):
    raw = np.zeros((4, 16, 16), dtype=np.uint16)
    tifffile.imwrite(tmp_path / "sq.ome.tif", raw, ome=True, metadata={
        "axes": "ZYX", "PhysicalSizeX": 0.0567, "PhysicalSizeY": 0.05670001, "PhysicalSizeZ": 0.3,
    })
    meta = get_metadata_from_tif(tmp_path / "sq.ome.tif")
    assert meta["pixel_size"] == pytest.approx(56.7) and "non_square_pixels" not in meta
