import numpy as np
from Batch import _interpolate_or_extrude_roi


def test_single_slice_extrusion():
    """1-slice ROI should extrude to all Z slices (legacy compatibility)."""
    mask = np.zeros((10, 50, 50), dtype=bool)
    mask[5, 20:30, 20:30] = True

    result = _interpolate_or_extrude_roi(mask)
    assert result.shape == (10, 50, 50)
    for z in range(10):
        assert np.array_equal(result[z], mask[5])


def test_multi_slice_interpolation_bounds():
    """Multi-slice ROI should interpolate between slices and be strictly zero outside."""
    mask = np.zeros((20, 50, 50), dtype=bool)
    # Small circle at z=5
    y, x = np.ogrid[:50, :50]
    mask[5] = (y - 25)**2 + (x - 25)**2 <= 5**2
    # Larger circle at z=15
    mask[15] = (y - 25)**2 + (x - 25)**2 <= 15**2

    result = _interpolate_or_extrude_roi(mask)

    # Outside [5, 15] must be strictly 0
    for z in range(0, 5):
        assert result[z].sum() == 0, f"Slice {z} should be empty"
    for z in range(16, 20):
        assert result[z].sum() == 0, f"Slice {z} should be empty"

    # Slices between 5 and 15 must be non-empty and monotonically increasing in area
    areas = [result[z].sum() for z in range(5, 16)]
    assert areas[0] > 0
    assert areas[-1] > areas[0]
    # Check that middle slice area is between slice 5 and slice 15
    assert areas[0] < areas[5] < areas[-1]


def test_empty_or_2d_graceful():
    """Empty or 2D arrays should pass through without error."""
    empty_3d = np.zeros((5, 10, 10), dtype=bool)
    assert np.array_equal(_interpolate_or_extrude_roi(empty_3d), empty_3d)

    arr_2d = np.ones((10, 10), dtype=bool)
    assert np.array_equal(_interpolate_or_extrude_roi(arr_2d), arr_2d)


# --- K3 (review 2026-09-23): several nuclei (IF) must stay separate -----------

import pytest
from scipy.ndimage import distance_transform_edt


def _legacy_interpolate(mask_3d):
    """Pre-K3 implementation, kept verbatim as the ExM (one nucleus) reference."""
    non_empty_z = np.where(mask_3d.any(axis=(1, 2)))[0]
    if len(non_empty_z) == 1:
        mask_2d = mask_3d.max(axis=0)
        return np.repeat(mask_2d[np.newaxis, :, :], mask_3d.shape[0], axis=0)
    out = np.zeros_like(mask_3d, dtype=mask_3d.dtype)
    for i in range(len(non_empty_z) - 1):
        z0, z1 = non_empty_z[i], non_empty_z[i + 1]
        m0, m1 = (mask_3d[z0] > 0), (mask_3d[z1] > 0)
        d0 = distance_transform_edt(m0) - distance_transform_edt(~m0)
        d1 = distance_transform_edt(m1) - distance_transform_edt(~m1)
        for z in range(z0, z1 + 1):
            w1 = (z - z0) / (z1 - z0)
            out[z] = ((1.0 - w1) * d0 + w1 * d1 > 0).astype(mask_3d.dtype)
    return out


def _ellipse(z, cy, cx, ry, rx):
    return np.array([[z, cy - ry, cx - rx], [z, cy - ry, cx + rx],
                     [z, cy + ry, cx + rx], [z, cy + ry, cx - rx]], dtype=float)


def _to_labels(shapes, shape=(20, 80, 80)):
    napari = pytest.importorskip("napari")
    layer = napari.layers.Shapes(shapes, shape_type="ellipse", ndim=3)
    return layer.to_labels(labels_shape=shape)


def test_exm_single_nucleus_matches_legacy_result():
    """(i) One nucleus drawn on 2–4 slices through napari: output identical to pre-K3."""
    cases = [
        [_ellipse(5, 40, 40, 8, 10), _ellipse(15, 40, 42, 20, 16)],
        [_ellipse(2, 40, 40, 5, 5), _ellipse(9, 38, 41, 18, 22), _ellipse(17, 42, 40, 7, 9)],
        [_ellipse(3, 30, 30, 6, 6), _ellipse(6, 32, 33, 12, 10),
         _ellipse(12, 35, 35, 15, 15), _ellipse(18, 36, 34, 4, 5)],
    ]
    for shapes in cases:
        labels = _to_labels(shapes)
        result = _interpolate_or_extrude_roi(labels)
        np.testing.assert_array_equal(result > 0, _legacy_interpolate(labels) > 0)
        assert set(np.unique(result)) == {0, 1}


def test_two_nuclei_each_on_two_slices_get_own_id_and_z_range():
    """(ii) Built through real napari Shapes.to_labels (shape IDs, not nucleus IDs)."""
    shapes = [
        _ellipse(3, 20, 20, 8, 8),   # A
        _ellipse(12, 60, 60, 8, 8),  # B
        _ellipse(8, 20, 20, 10, 10), # A
        _ellipse(17, 60, 60, 10, 10),# B
    ]
    result = _interpolate_or_extrude_roi(_to_labels(shapes))

    ids = sorted(set(np.unique(result)) - {0})
    assert ids == [1, 2]
    z_a = np.where((result == 1).any(axis=(1, 2)))[0]
    z_b = np.where((result == 2).any(axis=(1, 2)))[0]
    assert z_a.tolist() == list(range(3, 9))
    assert z_b.tolist() == list(range(12, 18))
    # Each ID stays in its own XY region.
    assert not result[:, 40:, 40:][result[:, 40:, 40:] > 0].__contains__(1)
    assert not result[:, :40, :40][result[:, :40, :40] > 0].__contains__(2)


def test_two_nuclei_each_on_one_different_slice_are_both_extruded():
    """(iii) The reported bug: A on z=5, B on z=20 left slices 6–19 empty."""
    labels = _to_labels([_ellipse(5, 20, 20, 8, 8), _ellipse(15, 60, 60, 8, 8)])
    result = _interpolate_or_extrude_roi(labels)

    assert sorted(set(np.unique(result)) - {0}) == [1, 2]
    for z in range(result.shape[0]):
        assert (result[z] == 1).sum() == (labels[5] > 0).sum()
        assert (result[z] == 2).sum() == (labels[15] > 0).sum()


def test_overlapping_shapes_on_same_slice_raise():
    """(iv) Two overlapping shapes on one slice cannot be separated."""
    labels = _to_labels([_ellipse(10, 40, 40, 8, 8), _ellipse(10, 44, 46, 8, 8)])
    with pytest.raises(ValueError, match="same Z-slice"):
        _interpolate_or_extrude_roi(labels)


def test_same_slice_shapes_joined_through_another_slice_raise():
    """Two separate shapes on z=10 both overlapping one shape on z=4 → ambiguous group."""
    labels = _to_labels([
        _ellipse(4, 40, 40, 6, 20),
        _ellipse(10, 40, 26, 6, 6),
        _ellipse(10, 40, 54, 6, 6),
    ])
    with pytest.raises(ValueError, match="same Z-slice"):
        _interpolate_or_extrude_roi(labels)


def test_group_overlap_guard_raises():
    """Overlap between interpolated groups is impossible with overlap-based grouping
    (SDF blend stays inside the union of a group's drawn shapes), so the guard is
    tested directly."""
    from Batch import _assert_roi_groups_disjoint

    a = np.zeros((4, 10, 10), dtype=bool)
    b = np.zeros((4, 10, 10), dtype=bool)
    a[1, 2:5, 2:5] = True
    b[1, 4:7, 4:7] = True
    with pytest.raises(ValueError, match="overlap"):
        _assert_roi_groups_disjoint([a, b])
    b[1, 4, 4] = False
    _assert_roi_groups_disjoint([a, b])


def test_binary_mask_input_still_supported():
    """Painted/boolean ROI input (one value) is grouped by connected components."""
    mask = np.zeros((10, 50, 50), dtype=bool)
    mask[2, 5:15, 5:15] = True
    mask[7, 5:15, 5:15] = True
    mask[4, 30:40, 30:40] = True
    result = _interpolate_or_extrude_roi(mask)
    assert sorted(set(np.unique(result)) - {0}) == [1, 2]
    assert np.where((result == 1).any(axis=(1, 2)))[0].tolist() == list(range(2, 8))
    assert (result == 2).any(axis=(1, 2)).all()
