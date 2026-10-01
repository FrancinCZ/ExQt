# K1 (review 2026-09-23): K_part must not depend on ROI shape or noise through a
# pseudo camera offset; the offset is an explicit setting, padding voxels from
# alignment never enter statistics, and Auto-ROI does not report a K.
import numpy as np
import pytest
import tifffile

from Batch import process_condensates

OFFSET = 100.0
SHAPE = (12, 96, 96)
CENTER = (6, 48, 48)
NUCLEUS_RADII = (5.0, 30.0, 30.0)
BLOB_CENTERS = [(48 + dy, 48 + dx) for dy, dx in ((0, 0), (10, 4), (-9, 7), (5, -11), (-6, -8), (12, -3))]


def _synthetic(sigma, seed=0):
    """Offset 100, nucleoplasm +100, condensates +500, gel outside nucleus +10 → true K = 5."""
    z, y, x = np.indices(SHAPE, dtype=float)
    nucleus = (((z - CENTER[0]) / NUCLEUS_RADII[0]) ** 2
               + ((y - CENTER[1]) / NUCLEUS_RADII[1]) ** 2
               + ((x - CENTER[2]) / NUCLEUS_RADII[2]) ** 2) <= 1.0
    mask = np.zeros(SHAPE, dtype=np.uint8)
    for cy, cx in BLOB_CENTERS:
        mask[5:8][((y[5:8] - cy) ** 2 + (x[5:8] - cx) ** 2) <= 9] = 1
    image = np.full(SHAPE, OFFSET + 10.0)
    image[nucleus] = OFFSET + 100.0
    image[mask > 0] = OFFSET + 500.0
    image += np.random.default_rng(seed).normal(0.0, sigma, SHAPE)
    return np.clip(np.rint(image), 0, 65535).astype(np.uint16), mask


def _disk(radius, cy=48, cx=48):
    y, x = np.indices(SHAPE[1:], dtype=float)
    return ((y - cy) ** 2 + (x - cx) ** 2) <= radius ** 2


def _roi_loose():
    # One slice, extruded: covers the nucleus up to its rim at the condensate Z range.
    roi = np.zeros(SHAPE, dtype=np.int64)
    roi[6] = _disk(29)
    return roi


def _roi_tight_two_slices():
    # Tight ROI inside the nucleus on two slices → interpolated between z=3 and z=9.
    roi = np.zeros(SHAPE, dtype=np.int64)
    roi[3] = _disk(20)
    roi[9] = _disk(20) * 2
    return roi


def _run(tmp_path, image, mask, roi=None, name="img", **kwargs):
    tif, msk = tmp_path / f"{name}.tif", tmp_path / f"{name}_Mask.tif"
    tifffile.imwrite(tif, image)
    tifffile.imwrite(msk, mask)
    return process_condensates(
        tif_path=tif, mask_path=msk, mode="3d", min_voxels=1,
        auto_roi=roi is None,
        request_roi_func=(lambda shape, is_3d: roi) if roi is not None else None,
        pixel_size_nm=100.0, z_step_nm=300.0, signal_channel=0, dapi_channel=0,
        **kwargs,
    )


@pytest.mark.parametrize("sigma", [5, 15, 30])
@pytest.mark.parametrize("roi_factory", [_roi_loose, _roi_tight_two_slices], ids=["loose", "tight"])
def test_manual_roi_k_is_independent_of_roi_shape_and_noise(tmp_path, sigma, roi_factory):
    image, mask = _synthetic(sigma)
    df = _run(tmp_path, image, mask, roi_factory(), detector_offset_adu=OFFSET)
    assert len(df) == len(BLOB_CENTERS)
    assert df["K_valid"].all()
    assert df["partition_coefficient"].median() == pytest.approx(5.0, abs=0.5)
    assert set(df["K_offset_method"]) == {"explicit_setting"}
    assert set(df["K_offset_adu"]) == {OFFSET}


def test_default_offset_is_zero_and_recorded(tmp_path):
    image, mask = _synthetic(5)
    df = _run(tmp_path, image, mask, _roi_tight_two_slices())
    assert set(df["K_offset_adu"]) == {0.0}
    # (600) / (200) without any offset subtraction.
    assert df["partition_coefficient"].median() == pytest.approx(3.0, abs=0.1)


def test_negative_or_non_finite_offset_is_rejected(tmp_path):
    image, mask = _synthetic(5)
    for bad in (-1.0, float("nan")):
        with pytest.raises(ValueError):
            _run(tmp_path, image, mask, _roi_tight_two_slices(), detector_offset_adu=bad)


def test_auto_roi_reports_no_k(tmp_path):
    image, mask = _synthetic(5)
    df = _run(tmp_path, image, mask, None, detector_offset_adu=OFFSET)
    assert not df.empty
    assert df["partition_coefficient"].isna().all()
    assert df["nucleoplasm_mean_intensity"].isna().all()
    assert (~df["K_valid"].astype(bool)).all()
    assert set(df["K_invalid_reason"]) == {"auto_roi_no_nucleoplasm"}


def test_alignment_padding_does_not_change_k(tmp_path):
    """Padded (expand_canvas) stack with its valid-voxel sidecar gives the same K."""
    image, mask = _synthetic(15)
    pad = ((0, 0), (9, 0), (0, 7))
    roi = np.ones(SHAPE, dtype=np.int64)            # whole original FOV
    roi_padded = np.ones(np.pad(roi, pad).shape, dtype=np.int64)  # ROI drawn over padding too

    plain = _run(tmp_path, image, mask, roi, name="plain", detector_offset_adu=OFFSET)

    padded_dir = tmp_path / "aligned"
    padded_dir.mkdir()
    valid = np.pad(np.ones(SHAPE, dtype=np.uint8), pad)
    tifffile.imwrite(padded_dir / "img_alignment_valid.tif", valid)
    padded = _run(padded_dir, np.pad(image, pad), np.pad(mask, pad), roi_padded,
                  detector_offset_adu=OFFSET)

    k_plain = plain["partition_coefficient"].median()
    assert padded["partition_coefficient"].median() == pytest.approx(k_plain, rel=0.02)

    # Sensitivity check: without the sidecar the zero padding would bias K.
    (padded_dir / "img_alignment_valid.tif").unlink()
    biased = _run(padded_dir, np.pad(image, pad), np.pad(mask, pad), roi_padded,
                  detector_offset_adu=OFFSET)
    assert biased["partition_coefficient"].median() > k_plain * 1.02


def test_valid_mask_shape_mismatch_is_an_error(tmp_path):
    image, mask = _synthetic(5)
    tifffile.imwrite(tmp_path / "img_alignment_valid.tif", np.ones((12, 10, 10), dtype=np.uint8))
    with pytest.raises(ValueError, match="valid-voxel"):
        _run(tmp_path, image, mask, _roi_tight_two_slices(), detector_offset_adu=OFFSET)


def test_aligner_writes_valid_voxel_sidecar(tmp_path):
    from scipy import ndimage as ndi
    from stack_aligner import AlignmentConfig, align_tiff_stack

    y, x = np.indices((80, 92), dtype=float)
    base = (3200 * np.exp(-((y - 20) ** 2 + (x - 28) ** 2) / 85)
            + 2700 * np.exp(-((y - 56) ** 2 + (x - 64) ** 2) / 115)
            + 1900 * np.exp(-((y - 40) ** 2 + (x - 48) ** 2) / 50))
    positions = np.array([[0, 0], [1, -1], [2, -2], [2, -3]], dtype=float)
    stack = np.stack([ndi.shift(base + 400, p, order=1, mode="constant", prefilter=False)
                      for p in positions]).astype(np.uint16)
    source = tmp_path / "in" / "s.tif"
    source.parent.mkdir()
    tifffile.imwrite(source, stack, ome=True, metadata={"axes": "ZYX"})

    result = align_tiff_stack(source, tmp_path / "out", config=AlignmentConfig(expand_canvas=True))
    valid_path = tmp_path / "out" / "s_alignment_valid.tif"
    assert valid_path.is_file()
    valid = tifffile.imread(valid_path)
    with tifffile.TiffFile(result.output) as tif:
        aligned = tif.asarray()
    assert valid.shape == aligned.shape
    assert set(np.unique(valid)) == {0, 1}
    # Padding is exactly zero in the aligned data and flagged invalid.
    assert np.all(aligned[valid == 0] <= aligned.max() * 0.5)
    assert (valid == 0).any() and (valid == 1).mean() > 0.8


def test_worker_passes_offset_and_records_partitioning_metadata(tmp_path):
    import json
    import os
    import pandas as pd

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from App import AnalysisWorker

    image, mask = _synthetic(5)
    tifffile.imwrite(tmp_path / "img.tif", image)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    worker = AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(tmp_path / "out"), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": True, "review_each_image": False,
        "show_napari": False, "generate_reports": False, "pixel_size_nm": 100.0, "z_step_nm": 300.0,
        "calibration_source": "gui", "signal_channel": 0, "dapi_channel": 0,
        "detector_offset_adu": OFFSET,
    })
    worker.run()
    out = tmp_path / "out"
    df = pd.read_csv(next(out.glob("*_Output_Batch_3d.csv")))
    meta = json.loads(next(out.glob("*_metadata.json")).read_text(encoding="utf-8"))
    assert set(df["K_offset_adu"]) == {OFFSET}
    assert meta["partitioning"]["detector_offset_adu"] == OFFSET
    assert meta["partitioning"]["K_offset_method"] == "explicit_setting"
    assert meta["partitioning"]["auto_roi"] is True
