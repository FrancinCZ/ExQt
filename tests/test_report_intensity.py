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
