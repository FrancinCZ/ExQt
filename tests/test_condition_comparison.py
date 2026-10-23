# Condition comparison (supervisor feedback 2026-10-17): runs grouped by the top-level folder of the
# merge root are drawn side by side on the same fixed axes as the standard plots. The nucleus, not the
# object, is the unit: objects of one nucleus are not independent measurements.
import json

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from defaults import DEFAULT_SETTINGS
from postprocessing import build_condition_comparison_figure, condition_label, nucleus_summary

MIN_SIZE, MAX_SIZE = 0.0001, 2.0


def _primary(run_id, k, n_objects=(3, 5)):
    rows = []
    for file_index, n in enumerate(n_objects):
        for i in range(n):
            rows.append(dict(run_id=run_id, filename=f"img{file_index}.tif", cell_id=1, object_id=i + 1,
                             volume_bio_um3=0.001 * (i + 1), partition_coefficient=k + 0.1 * i, K_valid=True,
                             A_object=0.3, primary_qc_valid=True))
    return pd.DataFrame(rows)


def test_condition_is_top_level_folder():
    assert condition_label("Control/rep1") == "Control"
    assert condition_label("Control") == "Control"
    assert condition_label("run_Output_Batch_3d") == "run_Output_Batch_3d"


def test_nucleus_summary_uses_one_row_per_nucleus():
    primary = pd.concat([_primary("Control/rep1", 2.0), _primary("Control/rep2", 2.0), _primary("Treated", 4.0)])
    summary = nucleus_summary(primary)
    # Same filename in two replicate runs is two different nuclei.
    assert len(summary) == 6
    assert set(summary["condition"]) == {"Control", "Treated"}
    assert summary["n_objects"].tolist().count(5) == 3


def test_invalid_k_is_not_used():
    primary = _primary("Control", 2.0)
    primary["K_valid"] = False
    summary = nucleus_summary(primary)
    assert summary["median_K"].isna().all()


def _panel(fig, letter):
    return next(ax for ax in fig.axes if ax.get_title().startswith(f"{letter})"))


def test_figure_has_fixed_axes_and_conditions_in_order():
    primary = pd.concat([_primary("B_cond", 4.0), _primary("A_cond", 2.0)])
    fig = build_condition_comparison_figure(primary, MIN_SIZE, MAX_SIZE)
    assert np.allclose(_panel(fig, "A").get_ylim(), (MIN_SIZE, MAX_SIZE))
    for letter in "BC":
        assert np.allclose(_panel(fig, letter).get_ylim(), (DEFAULT_SETTINGS["plot_k_min"], DEFAULT_SETTINGS["plot_k_max"]))
    for letter in "DE":
        assert np.allclose(_panel(fig, letter).get_ylim(), (0.0, 1.0))
    for letter in "CE":
        assert np.allclose(_panel(fig, letter).get_xlim(), (MIN_SIZE, MAX_SIZE))
    ticks = [t.get_text() for t in _panel(fig, "B").get_xticklabels()]
    assert ticks[0].startswith("A_cond") and ticks[1].startswith("B_cond")
    assert "2 nuclei" in ticks[0]


# Supervisor request 2026-10-23: Control / Test / Test2 in one figure, each of the standard graphs with one
# box per sample (like the size panel), and the scatter panels coloured by sample.
def test_six_standard_panels_with_one_box_per_condition():
    primary = pd.concat([_primary("Control", 2.0), _primary("Test", 4.0), _primary("Test2", 6.0)])
    fig = build_condition_comparison_figure(primary, MIN_SIZE, MAX_SIZE)
    titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    assert [t[:2] for t in titles] == ["A)", "B)", "C)", "D)", "E)", "F)"]
    for letter in "ABD":
        boxes = [p for p in _panel(fig, letter).patches if p.__class__.__name__ == "PathPatch"]
        assert len(boxes) == 3
    legend = _panel(fig, "C").get_legend()
    assert [t.get_text() for t in legend.get_texts()][:3] == ["Control", "Test", "Test2"]


def _write_run(folder, k):
    folder.mkdir(parents=True)
    df = _primary("x", k).drop(columns=["run_id", "primary_qc_valid"]).assign(
        mode="3d", is_3d=True, volume_px=600, touches_image_edge=False, touches_roi_edge=False)
    csv = folder / "run_Output_Batch_3d.csv"
    df.to_csv(csv, index=False)
    metadata = {"parameters": {"plot_min_size": MIN_SIZE, "plot_max_size": MAX_SIZE, "mode": "3d"}}
    (folder / "run_Output_Batch_3d_metadata.json").write_text(json.dumps(metadata), encoding="utf-8")


def test_merge_writes_condition_comparison(tmp_path):
    from postprocessing import merge_statistics_folder

    _write_run(tmp_path / "Control" / "rep1", 2.0)
    _write_run(tmp_path / "Treated" / "rep1", 4.0)
    result = merge_statistics_folder(tmp_path)
    assert result["comparison_plot"].is_file()
    nuclei = pd.read_excel(result["excel"], sheet_name="Per_Nucleus")
    assert set(nuclei["condition"]) == {"Control", "Treated"}
