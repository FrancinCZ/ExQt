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
        win.compare_runs()
    assert out.is_file()
    assert "IF: 5 nuclei" in info.call_args[0][2]
    win.close()
    app.processEvents()
