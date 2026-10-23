# Merge overview plot axes (feedback 2026-10-17): the Delta A panels printed tick values such as
# -0.05, 0.02, 0.09 because ticks were evenly spaced between data-driven limits. Ticks must be round
# numbers, include zero when the range spans it, and labels must say "1 nucleus", not "1 nuclei".
import numpy as np
import pandas as pd

from postprocessing import _count_text, _nice_axis, build_condition_comparison_figure


def test_nice_axis_gives_round_ticks_including_zero():
    low, high, ticks = _nice_axis(-0.05, 0.30)
    assert low <= -0.05 and high >= 0.30
    assert 0.0 in [round(t, 10) for t in ticks]
    steps = np.diff(ticks)
    assert np.allclose(steps, steps[0])
    assert all(abs(t / steps[0] - round(t / steps[0])) < 1e-9 for t in ticks)
    assert 3 <= len(ticks) <= 8
    assert ticks[0] == low and ticks[-1] == high


def test_nice_axis_unit_range_is_unchanged():
    low, high, ticks = _nice_axis(0.0, 1.0)
    assert (low, high) == (0.0, 1.0)
    assert np.allclose(ticks, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


def test_count_text_singular_and_plural():
    assert _count_text(1, "nucleus", "nuclei") == "1 nucleus"
    assert _count_text(3, "nucleus", "nuclei") == "3 nuclei"


def test_condition_ticks_use_singular_for_one_nucleus():
    primary = pd.DataFrame(dict(run_id="img001/Result", filename="a.tif", cell_id=1, object_id=[1, 2],
                                volume_bio_um3=[0.5, 0.6], partition_coefficient=[5.0, 6.0], K_valid=True,
                                A_object=0.4, primary_qc_valid=True))
    fig = build_condition_comparison_figure(primary, 0.2, 2.0)
    k_panel = next(ax for ax in fig.axes if ax.get_title().startswith("B)"))
    ticks = [t.get_text() for t in k_panel.get_xticklabels()]
    assert "1 nucleus /" in ticks[0]


def test_merge_plot_ticks_stay_inside_axis(tmp_path):
    from postprocessing import _generate_merge_plot

    run_stats = pd.DataFrame(dict(run_id=["a", "b"], median_delta_a=[0.2, 0.28], qc_acceptance_rate=[0.8, 0.7],
                                  primary_objects=[8, 7], size_eligible_objects=[10, 10]))
    primary = pd.DataFrame(dict(run_id=["a"] * 3 + ["b"] * 3, A_shell=0.6, A_middle=0.7, A_core=0.8,
                                Delta_A_core_shell=[-0.08, 0.1, 0.42, 0.0, 0.2, 0.3]))
    low, high, ticks = _nice_axis(-0.11, 0.45)
    # Feeding already-round limits back must not widen them (the plot reuses them for drawing).
    assert _nice_axis(low, high, max_ticks=len(ticks))[:2] == (low, high)
    assert _generate_merge_plot(tmp_path / "m.xlsx", run_stats, primary).is_file()
