# Compare runs (feedback 2026-10-20): the supervisor wants ExM and IF graphs side by side. Merge pools runs
# and requires one QC policy, so ExM (ExF 10) and IF (ExF 1) cannot be merged. Compare runs is a separate
# tool: each chosen folder is one group, runs inside a group must share a QC policy, and settings that
# differ between groups are listed on the figure and in the workbook instead of blocking the comparison.
import json

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest

from postprocessing import QCPolicyMismatchError, compare_run_groups


def _write_run(folder, k, expansion_factor, min_size, max_size, allow_z_edge=False, cells=1):
    folder.mkdir(parents=True)
    rows = [dict(filename="img.tif", cell_id=cell, object_id=i + 1, volume_bio_um3=min_size * (2 + i),
                 partition_coefficient=k + 0.1 * i, K_valid=True, A_object=0.4, mode="3d", is_3d=True,
                 volume_px=600, touches_image_edge=False, touches_roi_edge=False)
            for cell in range(1, cells + 1) for i in range(4)]
    pd.DataFrame(rows).to_csv(folder / "run_Output_Batch_3d.csv", index=False)
    metadata = {"parameters": {"plot_min_size": min_size, "plot_max_size": max_size, "mode": "3d",
                               "expansion_factor": expansion_factor, "min_voxels": 50,
                               "allow_z_edge": allow_z_edge}}
    (folder / "run_Output_Batch_3d_metadata.json").write_text(json.dumps(metadata), encoding="utf-8")


def _groups(root, exm_k=5.0):
    _write_run(root / "SON-ExM" / "img001" / "Result", exm_k, 10.0, 0.02, 2.0)
    _write_run(root / "SON-ExM" / "img002" / "Result", exm_k + 1, 10.0, 0.02, 2.0)
    _write_run(root / "IF" / "img001" / "Result", 16.0, 1.0, 0.9, 4.7, cells=5)
    return [("ExM", root / "SON-ExM"), ("IF", root / "IF")]


def test_compare_exm_and_if_side_by_side(tmp_path):
    result = compare_run_groups(_groups(tmp_path), tmp_path / "out" / "ExM_vs_IF.png")
    assert result["plot"].is_file() and result["excel"].is_file()
    nuclei = pd.read_excel(result["excel"], sheet_name="Per_Nucleus")
    assert nuclei.groupby("condition").size().to_dict() == {"ExM": 2, "IF": 5}
    groups = pd.read_excel(result["excel"], sheet_name="Per_Group").set_index("group")
    assert groups.loc["IF", "nuclei"] == 5 and groups.loc["ExM", "runs"] == 2
    assert groups.loc["IF", "median_K_of_nuclei"] == pytest.approx(16.15)


def test_differences_between_groups_are_reported_not_blocking(tmp_path):
    result = compare_run_groups(_groups(tmp_path), tmp_path / "cmp.png")
    assert any("expansion_factor" in line and "ExM=10" in line and "IF=1" in line for line in result["differences"])
    sheet = pd.read_excel(result["excel"], sheet_name="Settings_Differences")
    assert "analysis.expansion_factor" in set(sheet["parameter"])


def test_runs_inside_a_group_must_share_qc_policy(tmp_path):
    groups = _groups(tmp_path)
    _write_run(tmp_path / "SON-ExM" / "img003" / "Result", 5.0, 10.0, 0.02, 2.0, allow_z_edge=True)
    with pytest.raises(QCPolicyMismatchError) as error:
        compare_run_groups(groups, tmp_path / "cmp.png")
    assert "ExM" in str(error.value) and "allow_z_edge" in str(error.value)


def test_groups_keep_chosen_order_and_size_axis_covers_all(tmp_path):
    import postprocessing

    captured = {}
    original = postprocessing.build_condition_comparison_figure

    def spy(primary, min_size, max_size, **kwargs):
        captured.update(range=(min_size, max_size), **kwargs)
        return original(primary, min_size, max_size, **kwargs)

    postprocessing.build_condition_comparison_figure = spy
    try:
        compare_run_groups(list(reversed(_groups(tmp_path))), tmp_path / "cmp.png")
    finally:
        postprocessing.build_condition_comparison_figure = original
    assert np.allclose(captured["range"], (0.02, 4.7))
    assert captured["order"] == ["IF", "ExM"]


def test_group_labels_must_be_unique_and_folders_must_have_runs(tmp_path):
    groups = _groups(tmp_path)
    with pytest.raises(ValueError, match="unique"):
        compare_run_groups([("A", groups[0][1]), ("A", groups[1][1])], tmp_path / "cmp.png")
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="empty"):
        compare_run_groups([groups[0], ("Empty", tmp_path / "empty")], tmp_path / "cmp.png")


def test_merge_still_refuses_exm_and_if(tmp_path):
    from postprocessing import merge_statistics_folder

    _groups(tmp_path)
    with pytest.raises(QCPolicyMismatchError):
        merge_statistics_folder(tmp_path)


def test_compare_runs_dialog_and_action(tmp_path):
    import os
    from unittest.mock import patch
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication, QDialogButtonBox, QWidget
    from PySide6.QtCore import QSettings

    app = QApplication.instance() or QApplication([])
    groups = _groups(tmp_path)
    from App import CompareRunsDialog
    dialog = CompareRunsDialog()
    ok = dialog.buttons.button(QDialogButtonBox.Ok)
    dialog.add_folder(str(groups[0][1]))
    assert not ok.isEnabled()
    dialog.add_folder(str(groups[1][1]))
    assert ok.isEnabled()
    dialog.table.item(0, 0).setText("ExM")
    assert dialog.groups() == [("ExM", str(groups[0][1])), ("IF", str(groups[1][1]))]

    with patch("napari.Viewer") as viewer:
        viewer.return_value.window._qt_window = QWidget()
        from App import ExQt
        win = ExQt()
    win.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)
    assert win.compare_runs_action.text() == "Compare runs..."
    out = tmp_path / "ExM_vs_IF.png"
    with patch("App.CompareRunsDialog") as dialog_cls, \
            patch("App.QFileDialog.getSaveFileName", return_value=(str(out), "")), \
            patch("App.QMessageBox.information") as info:
        dialog_cls.return_value.exec.return_value = 1
        dialog_cls.return_value.groups.return_value = [("ExM", str(groups[0][1])), ("IF", str(groups[1][1]))]
        dialog_cls.return_value.selected_nuclei.return_value = None
        dialog_cls.return_value.per_nucleus.return_value = False
        win.compare_runs()
    assert out.is_file()
    assert "IF: 5 nuclei" in info.call_args[0][2]
    win.close()
    app.processEvents()


# Picking nuclei (request 2026-10-25): the supervisor wants 2 ExM and 2 IF nuclei compared graph by graph.
# Compare runs lists every nucleus of the chosen folders; only ticked nuclei are used, and optionally each
# nucleus is its own column (coloured by its group).
def _box_colors(ax):
    return [tuple(p.get_facecolor()) for p in ax.patches if p.__class__.__name__ == "PathPatch"]


def test_nuclei_are_listed_per_group(tmp_path):
    from postprocessing import list_compare_nuclei

    nuclei = list_compare_nuclei(_groups(tmp_path))
    assert [n["group"] for n in nuclei].count("ExM") == 2 and [n["group"] for n in nuclei].count("IF") == 5
    first = nuclei[0]
    assert first["label"] == "ExM img001 C1" and first["n_objects"] == 4
    assert len({n["key"] for n in nuclei}) == 7


def test_only_picked_nuclei_are_compared(tmp_path):
    from postprocessing import list_compare_nuclei

    groups = _groups(tmp_path)
    nuclei = list_compare_nuclei(groups)
    picked = [n["key"] for n in nuclei if n["label"] in ("ExM img001 C1", "IF img001 C2", "IF img001 C4")]
    result = compare_run_groups(groups, tmp_path / "pick.png", nuclei=picked)
    per_nucleus = pd.read_excel(result["excel"], sheet_name="Per_Nucleus")
    assert len(per_nucleus) == 3
    assert {row["group"]: row["nuclei"] for row in result["groups"]} == {"ExM": 1, "IF": 2}


def test_each_nucleus_as_its_own_column(tmp_path):
    import postprocessing
    from postprocessing import list_compare_nuclei

    groups = _groups(tmp_path)
    labels = ["ExM img001 C1", "ExM img002 C1", "IF img001 C2", "IF img001 C4"]
    picked = [n["key"] for n in list_compare_nuclei(groups) if n["label"] in labels]
    figures = []
    original = postprocessing.build_condition_comparison_figure

    def spy(*args, **kwargs):
        figures.append(original(*args, **kwargs))
        return figures[-1]

    postprocessing.build_condition_comparison_figure = spy
    try:
        compare_run_groups(groups, tmp_path / "each.png", nuclei=picked, per_nucleus=True)
    finally:
        postprocessing.build_condition_comparison_figure = original
    k_panel = next(ax for ax in figures[0].axes if ax.get_title().startswith("B)"))
    ticks = [t.get_text().split("\n")[0] for t in k_panel.get_xticklabels()]
    assert ticks == labels
    colors = _box_colors(k_panel)
    assert len(colors) == 4 and colors[0] == colors[1] != colors[2] == colors[3]


def test_unknown_nucleus_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="nucle"):
        compare_run_groups(_groups(tmp_path), tmp_path / "x.png", nuclei=["nope"])


def test_dialog_lists_nuclei_and_returns_ticked(tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    from App import CompareRunsDialog
    groups = _groups(tmp_path)
    dialog = CompareRunsDialog()
    dialog.add_folder(str(groups[0][1]), "ExM")
    dialog.add_folder(str(groups[1][1]), "IF")
    assert dialog.nuclei_list.count() == 7
    assert dialog.selected_nuclei() is None          # all ticked = no filter
    for row in range(dialog.nuclei_list.count()):
        item = dialog.nuclei_list.item(row)
        if not item.text().startswith(("ExM img001 C1", "IF img001 C2")):
            item.setCheckState(Qt.Unchecked)
    assert len(dialog.selected_nuclei()) == 2
    dialog.per_nucleus_check.setChecked(True)
    assert dialog.per_nucleus()
    dialog.close()
    app.processEvents()


def test_nucleus_labels_stay_unique_when_runs_share_a_folder(tmp_path):
    # Argon C (2026-10-26): two single-file runs under one date folder were both called "ExM 2026.08.06 C1",
    # so "each nucleus as its own column" pooled two nuclei into one column.
    from postprocessing import list_compare_nuclei

    _write_run(tmp_path / "ExM" / "2026.08.06" / "img001" / "Result", 5.0, 10.0, 0.02, 2.0)
    _write_run(tmp_path / "ExM" / "2026.08.06" / "img002" / "Result", 6.0, 10.0, 0.02, 2.0)
    _write_run(tmp_path / "IF" / "img001" / "Result", 16.0, 1.0, 0.9, 4.7, cells=2)
    groups = [("ExM", tmp_path / "ExM"), ("IF", tmp_path / "IF")]

    nuclei = list_compare_nuclei(groups)
    labels = [n["label"] for n in nuclei]
    assert len(set(labels)) == len(labels) == 4
    # Labels that did not collide keep their short form.
    assert "IF img001 C1" in labels and "IF img001 C2" in labels

    import postprocessing

    captured = {}
    original = postprocessing.build_condition_comparison_figure

    def spy(primary, min_size, max_size, **kwargs):
        captured["columns"] = primary.groupby("condition").size().to_dict()
        return original(primary, min_size, max_size, **kwargs)

    postprocessing.build_condition_comparison_figure = spy
    try:
        compare_run_groups(groups, tmp_path / "each.png", per_nucleus=True)
    finally:
        postprocessing.build_condition_comparison_figure = original
    assert len(captured["columns"]) == 4 and set(captured["columns"].values()) == {4}
