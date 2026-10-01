# The left control panel must not push the Napari viewer aside: long status messages wrap,
# and the split between panel and viewer is adjustable and remembered.
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture
def window(tmp_path):
    from unittest.mock import patch
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QWidget

    app = QApplication.instance() or QApplication([])
    with patch("napari.Viewer") as viewer:
        viewer.return_value.window._qt_window = QWidget()
        from App import ExQt
        win = ExQt()
    win.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.IniFormat)  # never the real registry
    yield win
    win.close()
    app.processEvents()


def test_long_status_text_does_not_widen_the_left_panel(window):
    panel = window.main_splitter.widget(0)
    before = panel.minimumSizeHint().width()
    window.status_label.setText(
        "Draw one shape per nucleus on 1 or more Z-slices (overlapping shapes on different slices = same "
        "nucleus, interpolated; outside its slices = 0; one slice = whole Z). Shapes must not touch on the "
        "same slice. Click 'Confirm ROI'..." * 2
    )
    assert window.status_label.wordWrap()
    assert panel.minimumSizeHint().width() <= before + 5


def test_panel_and_viewer_are_in_an_adjustable_splitter(window):
    from PySide6.QtWidgets import QSplitter

    splitter = window.main_splitter
    assert isinstance(splitter, QSplitter) and window.centralWidget() is splitter
    assert splitter.count() == 2 and not splitter.childrenCollapsible()


def test_splitter_position_is_remembered(window):
    window.main_splitter.setSizes([300, 900])
    window.close()
    assert window.settings.value("main_splitter_state") is not None
