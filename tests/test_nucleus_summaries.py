# Unit = nucleus (ROI) in the PDF, Merge and Compare runs (supervisor feedback 2026-11-25):
# mean ± SD for means, median (IQR) for medians, class shares averaged over nuclei.
import numpy as np
import pandas as pd
import pytest

from nucleus_stats import (
    ALL_NUCLEI, class_totals_figure, nucleus_columns_table, per_nucleus_table, summary_tables,
)


def _objects():
    rows = []
    # Nucleus 1: 2 objects (Elongated/Core, Globular/Core); nucleus 2: 4 Globular/No gradient + 1 excluded;
    # nucleus 3: only an excluded object.
    for cid, cls, k, vol, primary in [
        (1, "Elongated (High FA) / Core-Enriched", 10.0, 0.1, True),
        (1, "Globular (Low FA) / Core-Enriched", 20.0, 0.3, True),
        (2, "Globular (Low FA) / No Significant Gradient", 5.0, 0.2, True),
        (2, "Globular (Low FA) / No Significant Gradient", 5.0, 0.2, True),
        (2, "Globular (Low FA) / No Significant Gradient", 7.0, 0.2, True),
        (2, "Globular (Low FA) / No Significant Gradient", 7.0, 0.2, True),
        (2, "Globular (Low FA) / Core-Enriched", 99.0, 9.0, False),
        (3, "Unclassified (QC)", 1.0, 0.1, False),
    ]:
        rows.append({"filename": "img.tif", "cell_id": cid, "condensate_class": cls, "partition_coefficient": k,
                     "K_valid": True, "volume_bio_um3": vol, "mean_intensity": 100.0 * k,
                     "nucleoplasm_mean_intensity": 100.0, "A_object": 0.5, "primary_qc_valid": primary})
    return pd.DataFrame(rows)


def test_per_nucleus_table_counts_only_primary_and_keeps_empty_nuclei():
    nuclei = per_nucleus_table(_objects()).set_index("cell_id")
    assert nuclei["n_condensates"].to_dict() == {1: 2, 2: 4, 3: 0}
    assert nuclei.loc[1, "mean_K"] == pytest.approx(15.0)
    assert nuclei.loc[2, "median_volume"] == pytest.approx(0.2)
    assert nuclei.loc[1, "Total Elongated (High FA)"] == pytest.approx(50.0)
    assert nuclei.loc[2, "Total No Significant Gradient"] == pytest.approx(100.0)
    assert np.isnan(nuclei.loc[3, "Total Globular (Low FA)"])


def test_summary_is_mean_sd_over_nuclei_and_class_shares_weigh_nuclei_equally():
    table, numbers = summary_tables(per_nucleus_table(_objects()))
    column = table[ALL_NUCLEI]
    assert column["Nuclei (ROIs)"] == "3"
    assert column["Condensates per nucleus (ROI)"] == "2.0 ± 2.0"  # 2, 4, 0
    # Nucleus means 15 and 6; the nucleus without objects has no K.
    assert column["Mean partitioning (K_part)"] == "10.50 ± 6.36"
    # 50 % and 0 % Elongated -> 25 %, although pooled objects would give 1/6.
    assert column["Total Elongated (High FA)"] == "25.0 ± 35.4%"
    assert column["Median particle volume [µm³]"].startswith("0.2000 (")
    row = numbers[(numbers["metric"] == "Mean partitioning (K_part)")].iloc[0]
    assert row["n"] == 2 and row["summary"] == "mean ± SD"


def test_groups_side_by_side_in_given_order():
    objects = _objects()
    objects["group"] = np.where(objects["cell_id"] == 1, "HPV16", "pBabe")
    nuclei = per_nucleus_table(objects, group_cols=["group"])
    table, _ = summary_tables(nuclei, "group", ["pBabe", "HPV16"])
    assert list(table.columns) == ["pBabe", "HPV16"]
    assert table.loc["Nuclei (ROIs)", "pBabe"] == "2"
    assert table.loc["Condensates per nucleus (ROI)", "HPV16"] == "2.0"
    fig = class_totals_figure(nuclei, "group", ["pBabe", "HPV16"], {"pBabe": "#1f77b4", "HPV16": "#ff7f0e"})
    assert len(fig.axes) == 2


def test_nucleus_columns_table_has_one_column_per_nucleus_and_all():
    table = nucleus_columns_table(per_nucleus_table(_objects()))
    assert list(table.columns) == ["img / Cell 1", "img / Cell 2", "img / Cell 3", ALL_NUCLEI]
    assert table.loc["Condensates per nucleus (ROI)", "img / Cell 2"] == "4"
    assert table.loc["Total Elongated (High FA)", "img / Cell 1"] == "50.0%"


def test_pdf_has_per_nucleus_page(tmp_path):
    from matplotlib.backends.backend_pdf import PdfPages  # noqa: F401
    from partitioning_plots import export_partitioning_analysis

    objects = _objects()
    objects["object_id"] = range(1, len(objects) + 1)
    export_partitioning_analysis(objects, tmp_path, file_stem="run", min_size=0.01, max_size=10.0)
    pdf = (tmp_path / "run_ExQt_Report.pdf").read_bytes()
    assert pdf.count(b"/Type /Page\n") + pdf.count(b"/Type /Page ") + pdf.count(b"/Type/Page") >= 3


def test_merge_and_compare_write_nucleus_outputs(tmp_path):
    import json

    from postprocessing import compare_run_groups, merge_statistics_folder

    def write_run(folder, k):
        folder.mkdir(parents=True)
        objects = _objects()
        objects["partition_coefficient"] *= k
        objects["object_id"] = range(1, len(objects) + 1)
        for column in ("A_object_valid", "A_shell_valid", "A_middle_valid", "A_core_valid", "mode_a_primary_include"):
            objects[column] = True
        objects["mode_a_primary_include"] = objects.pop("primary_qc_valid")
        objects["A_shell"] = objects["A_middle"] = objects["A_core"] = 0.5
        csv = folder / "run_Output_Batch_3d.csv"
        objects.to_csv(csv, index=False)
        (folder / "run_Output_Batch_3d_metadata.json").write_text(json.dumps(
            {"parameters": {"plot_min_size": 0.01, "plot_max_size": 10.0, "expansion_factor": 1.0}}))
        return csv

    write_run(tmp_path / "pBabe" / "r1", 1.0)
    write_run(tmp_path / "HPV16" / "r1", 2.0)
    merged = merge_statistics_folder(tmp_path, tmp_path / "Merged_Stats.xlsx")
    sheets = pd.read_excel(merged["excel"], sheet_name=None, index_col=None)
    summary = sheets["Summary_Per_Nucleus"].set_index("Metric")
    assert list(summary.columns) == ["HPV16", "pBabe"]
    assert str(summary.loc["Nuclei (ROIs)", "pBabe"]) == "3"
    for key in ("summary_plot", "nuclei_plot", "class_totals_plot"):
        assert merged[key] is not None and merged[key].stat().st_size > 0

    result = compare_run_groups([("pBabe", tmp_path / "pBabe"), ("HPV16", tmp_path / "HPV16")],
                                tmp_path / "cmp" / "Compare")
    summary = pd.read_excel(result["excel"], sheet_name="Summary_Per_Nucleus").set_index("Metric")
    assert list(summary.columns) == ["pBabe", "HPV16"]
    assert summary.loc["Mean partitioning (K_part)", "HPV16"] == "21.00 ± 12.73"
    assert result["summary_plot"].stat().st_size > 0 and result["class_totals_plot"].stat().st_size > 0
