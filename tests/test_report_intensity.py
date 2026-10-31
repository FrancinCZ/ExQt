# PDF summary table (question 2026-10-21): the table of QC-passed particles should also give the mean
# intensity. It is raw ADU (comparable only between images acquired with the same settings), so the
# nucleoplasm level it is divided by for K_part is shown next to it.
import numpy as np
import pandas as pd
import pytest

from partitioning_plots import _summary_table_rows


def _frame():
    return pd.DataFrame(dict(
        primary_qc_valid=True, partition_coefficient=[4.0, 6.0, 8.0], K_valid=True,
        volume_bio_um3=[0.1, 0.2, 0.3], mean_intensity=[400.0, 600.0, 1100.0],
        nucleoplasm_mean_intensity=[100.0, 100.0, 120.0],
    ))


def _row(rows, label):
    return dict(rows)[label]


def test_table_has_mean_intensity_and_nucleoplasm():
    rows = _summary_table_rows(_frame(), "s")
    assert _row(rows, "Median Particle Mean Intensity") == "600.0 ADU"
    assert _row(rows, "Mean Particle Mean Intensity") == "700.0 ADU"
    assert _row(rows, "Median Nucleoplasm Intensity") == "100.0 ADU"
    assert _row(rows, "Median Partitioning (K_part)") == "6.00"


def test_missing_intensity_is_na():
    rows = _summary_table_rows(_frame().drop(columns=["mean_intensity", "nucleoplasm_mean_intensity"]), "s")
    assert _row(rows, "Median Particle Mean Intensity") == "N/A"
    assert _row(rows, "Median Nucleoplasm Intensity") == "N/A"


def test_summary_json_has_median_intensity(tmp_path):
    from partitioning_plots import export_partitioning_analysis

    summary = export_partitioning_analysis(_frame(), tmp_path, "s", export_pngs=False)
    assert summary["median_mean_intensity_adu"] == pytest.approx(600.0)
    assert (tmp_path / "s_ExQt_Report.pdf").is_file()


def test_class_rows_have_one_fixed_order_with_zeros():
    # Supervisor feedback 2026-10-31: class rows sorted by count changed order between PDFs.
    # The six main classes now always appear in one fixed order (0 included), others after them.
    classes = (["Globular (Low FA) / Shell-Enriched"] * 3 + ["Elongated (High FA) / No Significant Gradient"] * 2
               + ["Unclassified", "Elongated (High FA) / Core-Enriched"])
    frame = pd.DataFrame(dict(primary_qc_valid=True, volume_bio_um3=0.1, particle_class=classes))
    class_rows = [(label, value) for label, value in _summary_table_rows(frame, "s") if label.startswith("Class: ")]
    assert class_rows == [
        ("Class: Elongated (High FA) / Core-Enriched", "1 (14.3%)"),
        ("Class: Globular (Low FA) / Core-Enriched", "0 (0.0%)"),
        ("Class: Elongated (High FA) / No Significant Gradient", "2 (28.6%)"),
        ("Class: Globular (Low FA) / No Significant Gradient", "0 (0.0%)"),
        ("Class: Elongated (High FA) / Shell-Enriched", "0 (0.0%)"),
        ("Class: Globular (Low FA) / Shell-Enriched", "3 (42.9%)"),
        ("Class: Unclassified", "1 (14.3%)"),
    ]


def test_no_class_rows_without_classes():
    assert not [row for row in _summary_table_rows(_frame(), "s") if row[0].startswith("Class: ")]


def test_unclassified_only_run_has_no_zero_rows():
    frame = pd.DataFrame(dict(primary_qc_valid=True, volume_bio_um3=0.1, particle_class=["Unclassified"] * 2))
    class_rows = [row for row in _summary_table_rows(frame, "s") if row[0].startswith("Class: ")]
    assert class_rows == [["Class: Unclassified", "2 (100.0%)"]]
