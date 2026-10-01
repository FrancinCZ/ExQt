# K6 (review 2026-09-23): the K_part/classification report must only contain
# QC-passed objects, and the core/shell gradient class must not react to noise.
import json
import os

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from rezim_a_metrics import GRADIENT_ALPHA, classify_particle, compute_core_shell_metrics

SAMPLING = (100.0, 100.0, 100.0)


def _sphere(radius=8, pad=2):
    size = 2 * (radius + pad) + 1
    z, y, x = np.indices((size, size, size), dtype=float) - (radius + pad)
    return (z ** 2 + y ** 2 + x ** 2) <= radius ** 2


def _class_of(mask, intensity, **kwargs):
    return compute_core_shell_metrics(
        mask, sampling=SAMPLING, intensity_image=intensity, min_core_voxels=20, **kwargs
    )


def test_homogeneous_noisy_object_has_no_significant_gradient():
    """Acceptance: homogeneous object + noise → 'No Significant Gradient' in ≥ 90 %."""
    mask = _sphere()
    rng = np.random.default_rng(0)
    classes = []
    for _ in range(200):
        intensity = 500.0 + rng.normal(0.0, 15.0, mask.shape)
        classes.append(_class_of(mask, intensity)["condensate_class"])
    no_gradient = np.mean([c.endswith("No Significant Gradient") for c in classes])
    assert no_gradient >= 0.90
    # With alpha = 0.01 the expected false-positive rate is ~1 %.
    assert no_gradient >= 0.95


def test_real_core_enrichment_is_detected():
    mask = _sphere()
    z, y, x = np.indices(mask.shape, dtype=float) - 10
    r = np.sqrt(z ** 2 + y ** 2 + x ** 2)
    intensity = 500.0 + 60.0 * (1.0 - r / 8.0) + np.random.default_rng(1).normal(0.0, 15.0, mask.shape)
    metrics = _class_of(mask, intensity)
    assert metrics["condensate_class"].endswith("Core-Enriched")
    assert metrics["gradient_p_core_shell"] < GRADIENT_ALPHA
    assert metrics["gradient_test"] == "welch_t_two_sided"


def test_noise_free_constant_object_is_not_a_gradient():
    mask = _sphere()
    metrics = _class_of(mask, np.full(mask.shape, 500.0))
    assert metrics["condensate_class"].endswith("No Significant Gradient")


def test_qc_failed_object_is_unclassified():
    mask = _sphere()
    metrics = _class_of(mask, np.full(mask.shape, 500.0), primary_include=False,
                        primary_exclusion_reason="touches_image_edge")
    assert metrics["condensate_class"] == "Unclassified (QC)"
    assert classify_particle(0.1, 50.0, 5000, gradient_significant=True, qc_passed=False) == "Unclassified (QC)"


def _report_frame():
    """Two otherwise identical objects; the second failed Mode A QC."""
    base = {
        "filename": "a.tif", "cell_id": 1, "volume_bio_um3": 0.5, "volume_px": 500,
        "mean_intensity": 600.0, "partition_coefficient": 5.0, "K_valid": True,
        "A_object": 0.3, "A_shell": 0.3, "A_middle": 0.3, "A_core": 0.3,
        "A_object_valid": True, "A_shell_valid": True, "A_middle_valid": True, "A_core_valid": True,
        "condensate_class": "Globular (Low FA) / Core-Enriched", "mode_a_primary_include": True,
    }
    excluded = dict(base, object_id=2, mode_a_primary_include=False,
                    condensate_class="Unclassified (QC)", partition_coefficient=50.0)
    k_invalid = dict(base, object_id=3, K_valid=False, partition_coefficient=99.0)
    return pd.DataFrame([dict(base, object_id=1), excluded, k_invalid])


def test_report_contains_only_qc_passed_objects(tmp_path, monkeypatch):
    import partitioning_plots
    from postprocessing import _prepare_reporting_frames

    seen = {}
    real_pdf = partitioning_plots.export_unified_pdf_report

    def spy(df, *args, **kwargs):
        seen["pdf_ids"] = sorted(df["object_id"])
        return real_pdf(df, *args, **kwargs)

    monkeypatch.setattr(partitioning_plots, "export_unified_pdf_report", spy)
    _, _, _, report_df, _ = _prepare_reporting_frames(_report_frame(), 0.0001, 2.0)
    summary = partitioning_plots.export_partitioning_analysis(report_df, tmp_path, "s", export_pngs=False)

    assert seen["pdf_ids"] == [1, 3]                    # QC-failed object 2 never reaches the PDF
    assert summary["total_particles"] == 2
    assert "Unclassified (QC)" not in summary["class_percentages"]
    assert summary["median_partition_coefficient"] == pytest.approx(5.0)  # K_valid=False object 3 ignored
    assert "PSF" in summary["interpretation_warning"]
    saved = json.loads((tmp_path / "s_partitioning_summary.json").read_text(encoding="utf-8"))
    assert saved["total_particles"] == 2


def test_raw_frame_without_qc_columns_is_rejected(tmp_path):
    from partitioning_plots import export_partitioning_analysis

    with pytest.raises(ValueError, match="primary_qc_valid"):
        export_partitioning_analysis(_report_frame(), tmp_path, "s", export_pngs=False, export_pdf=False)


def test_worker_partitioning_report_uses_qc_frame(tmp_path):
    """Small objects fail Mode A QC (core too small) → they must not appear in the summary."""
    import tifffile
    from App import AnalysisWorker

    raw = np.full((8, 40, 40), 200, dtype=np.uint16)
    mask = np.zeros((8, 40, 40), dtype=np.uint8)
    for cy, cx in ((10, 10), (10, 28), (28, 18)):
        mask[3:5, cy - 2:cy + 2, cx - 2:cx + 2] = 1
        raw[3:5, cy - 2:cy + 2, cx - 2:cx + 2] = 900
    tifffile.imwrite(tmp_path / "img.tif", raw)
    tifffile.imwrite(tmp_path / "img_Mask.tif", mask)
    out = tmp_path / "out"
    worker = AnalysisWorker({
        "input_folder": str(tmp_path), "output_folder": str(out), "mode": "3d",
        "expansion_factor": 1.0, "min_voxels": 1, "auto_roi": True, "review_each_image": False,
        "show_napari": False, "generate_reports": True, "report_excel": False,
        "report_primary_csv": False, "report_excluded_csv": False, "report_raw_audit_csv": False,
        "report_standard_plots": False, "report_mode_a_plots": False, "report_partitioning_plots": True,
        "plot_min_size": 0.0, "plot_max_size": 1e6,
        "pixel_size_nm": 100.0, "z_step_nm": 300.0, "calibration_source": "gui",
        "signal_channel": 0, "dapi_channel": 0, "mode_a_enabled": True, "mode_a_min_core_voxels": 20,
    })
    messages = []
    worker.progress.connect(messages.append)
    worker.run()
    raw_df = pd.read_csv(next(out.glob("*_Output_Batch_3d.csv")))
    assert len(raw_df) == 3 and not raw_df["mode_a_primary_include"].any()
    summary = json.loads(next(out.glob("*_partitioning_summary.json")).read_text(encoding="utf-8"))
    assert summary["total_particles"] == 0, messages
