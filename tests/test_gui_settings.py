# One "Settings" entry opens the settings window directly; dark mode lives there and is
# remembered; the control panel stays compact.
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
    app.setStyleSheet("")
    app.processEvents()


def test_settings_is_a_single_menu_bar_action(window):
    from unittest.mock import patch
    import App

    actions = {a.text(): a for a in window.menuBar().actions()}
    assert "Settings" in actions and actions["Settings"].menu() is None   # no sub-menu
    with patch.object(App.AdvancedSettingsDialog, "exec") as exec_:
        actions["Settings"].trigger()
    assert exec_.called


def test_dark_mode_is_set_in_the_settings_window_and_remembered(window, tmp_path):
    from PySide6.QtWidgets import QApplication
    from App import AdvancedSettingsDialog

    dialog = AdvancedSettingsDialog(window)
    dialog.settings = window.settings
    dialog.dark_mode_check.setChecked(True)
    dialog.accept()
    assert window.settings.value("dark_mode") in (True, "true")
    dark = QApplication.instance().styleSheet()
    assert dark

    QApplication.instance().setStyleSheet("")
    window.load_settings()                      # next start: the saved choice is applied again
    assert QApplication.instance().styleSheet() == dark


def test_control_panel_is_compact(window):
    # Offscreen font metrics are about twice the real ones; before this change the minimum was 848.
    assert window.main_splitter.widget(0).minimumSizeHint().width() <= 650


def test_size_unit_follows_process_mode(window):
    window.mode_combo.setCurrentText("2d")
    assert window.plot_min_size_spin.suffix() == " µm²"
    window.mode_combo.setCurrentText("3d")
    assert window.plot_max_size_spin.suffix() == " µm³"
