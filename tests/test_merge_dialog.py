# Merge runs window (2026-11-07): one folder, a table of the runs found in it with their condition,
# nuclei, objects and QC differences; unticked runs are left out of the merge.
import json
import os

import pandas as pd

MIN_SIZE, MAX_SIZE = 0.0001, 2.0


def _write_run(folder, k, expansion=8.0):
    folder.mkdir(parents=True)
    rows = [{"filename": f"f{i}.tif", "cell_id": 1, "object_id": j, "volume_bio_um3": 0.1 + 0.01 * j,
             "volume_px": 600, "partition_coefficient": k, "K_valid": True, "mode": "3d", "is_3d": True,
             "touches_image_edge": False, "touches_roi_edge": False}
            for i in range(2) for j in range(3)]
    pd.DataFrame(rows).to_csv(folder / "run_Output_Batch_3d.csv", index=False)
    metadata = {"parameters": {"plot_min_size": MIN_SIZE, "plot_max_size": MAX_SIZE, "mode": "3d",
                               "expansion_factor": expansion}}
    (folder / "run_Output_Batch_3d_metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    return folder / "run_Output_Batch_3d.csv"


def _tree(tmp_path):
    a = _write_run(tmp_path / "Control" / "rep1", 2.0)
    b = _write_run(tmp_path / "Control" / "rep2", 2.5)
    c = _write_run(tmp_path / "Treated" / "rep1", 4.0, expansion=4.0)
    return a, b, c


def test_list_merge_runs_marks_the_odd_run(tmp_path):
    from postprocessing import list_merge_runs
    _tree(tmp_path)
    runs = {run["run_id"]: run for run in list_merge_runs(tmp_path)}
    assert set(runs) == {"Control/rep1", "Control/rep2", "Treated/rep1"}
    assert runs["Control/rep1"]["condition"] == "Control"
    assert runs["Control/rep1"]["nuclei"] == 2 and runs["Control/rep1"]["objects"] == 6
    assert runs["Control/rep1"]["differences"] == []
    assert runs["Treated/rep1"]["differences"]


def test_merge_only_selected_runs(tmp_path):
    from postprocessing import merge_statistics_folder
    a, b, _ = _tree(tmp_path)
    result = merge_statistics_folder(tmp_path, csv_paths=[a, b])
    assert result["run_count"] == 2
    assert set(pd.read_excel(result["excel"], sheet_name="Per_Run_Stats")["run_id"]) == {"Control/rep1", "Control/rep2"}


def test_merge_dialog_and_action(tmp_path):
    from unittest.mock import patch
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import QSettings, Qt
    from PySide6.QtWidgets import QApplication, QDialogButtonBox, QWidget

    app = QApplication.instance() or QApplication([])
    a, b, c = _tree(tmp_path)
    from App import MergeRunsDialog
    dialog = MergeRunsDialog()
    ok = dialog.buttons.button(QDialogButtonBox.Ok)
    assert not ok.isEnabled()
    dialog.set_folder(str(tmp_path))
    assert dialog.table.rowCount() == 3
    assert ok.text() == "Merge" and ok.isEnabled()
    assert sorted(dialog.selected_csvs()) == sorted([str(a), str(b)])
    assert dialog.table.item(2, 4).text().startswith("Differs:")
    for row in range(dialog.table.rowCount()):
        dialog.table.item(row, 0).setCheckState(Qt.Unchecked)
    assert not ok.isEnabled()
    assert not dialog.include_raw()

    with patch("napari.Viewer") as viewer:
        viewer.return_value.window._qt_window = QWidget()
        from App import ExQt
        win = ExQt()
    win.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)
    out = tmp_path / "Merged_Stats.xlsx"
    with patch("App.MergeRunsDialog") as dialog_cls, \
            patch("App.QFileDialog.getSaveFileName", return_value=(str(out), "")), \
            patch("App.QMessageBox.information") as info:
        dialog_cls.return_value.exec.return_value = 1
        dialog_cls.return_value.root = str(tmp_path)
        dialog_cls.return_value.selected_csvs.return_value = [str(a), str(b)]
        dialog_cls.return_value.include_raw.return_value = False
        win.merge_existing_runs()
    assert out.is_file()
    assert "Merged 2 compatible runs" in info.call_args[0][2]
    win.close()
