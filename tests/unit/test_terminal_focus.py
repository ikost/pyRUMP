"""OS-specific, best-effort focus handoff -- mocked so no real osascript,
xdotool, or user32 call ever runs in CI."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from pyrump.shell import terminal_focus


def _fake_run(stdout: str):
    def run(*args, **kwargs):
        return SimpleNamespace(stdout=stdout, returncode=0)

    return run


def test_capture_returns_none_on_an_unhandled_platform(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "freebsd13")
    assert terminal_focus.capture() is None


def test_restore_with_no_token_calls_nothing(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "darwin")

    def boom(*args, **kwargs):
        raise AssertionError("restore(None) should not shell out")

    monkeypatch.setattr(terminal_focus.subprocess, "run", boom)
    terminal_focus.restore(None)  # no exception means it never called run()


def test_macos_capture_reads_the_frontmost_process_name(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "darwin")
    monkeypatch.setattr(terminal_focus.subprocess, "run", _fake_run("iTerm2\n"))
    assert terminal_focus.capture() == "iTerm2"


def test_macos_restore_activates_the_captured_app(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "darwin")
    calls = []
    monkeypatch.setattr(
        terminal_focus.subprocess, "run",
        lambda args, **kwargs: calls.append(args) or SimpleNamespace(returncode=0),
    )
    terminal_focus.restore("iTerm2")
    assert calls and "iTerm2" in calls[0][-1]


def test_macos_capture_swallows_errors(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "darwin")

    def raises(*args, **kwargs):
        raise OSError("osascript not found")

    monkeypatch.setattr(terminal_focus.subprocess, "run", raises)
    assert terminal_focus.capture() is None


def test_macos_restore_swallows_errors(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "darwin")

    def raises(*args, **kwargs):
        raise OSError("osascript not found")

    monkeypatch.setattr(terminal_focus.subprocess, "run", raises)
    terminal_focus.restore("iTerm2")  # must not raise


def test_windows_capture_and_restore_use_user32(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "win32")
    calls = []
    fake_user32 = SimpleNamespace(
        GetForegroundWindow=lambda: 4242,
        SetForegroundWindow=lambda hwnd: calls.append(hwnd),
    )
    fake_windll = SimpleNamespace(user32=fake_user32)
    monkeypatch.setattr("ctypes.windll", fake_windll, raising=False)

    token = terminal_focus.capture()
    assert token == 4242
    terminal_focus.restore(token)
    assert calls == [4242]


def test_windows_capture_swallows_errors_when_windll_is_absent(monkeypatch):
    """``ctypes.windll`` doesn't exist off Windows -- capture must not raise."""
    monkeypatch.setattr(terminal_focus.sys, "platform", "win32")
    monkeypatch.delattr("ctypes.windll", raising=False)
    assert terminal_focus.capture() is None


def test_linux_capture_and_restore_use_xdotool(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "linux")
    monkeypatch.setattr(terminal_focus.subprocess, "run", _fake_run("0x1234\n"))
    assert terminal_focus.capture() == "0x1234"

    calls = []
    monkeypatch.setattr(
        terminal_focus.subprocess, "run",
        lambda args, **kwargs: calls.append(args) or SimpleNamespace(returncode=0),
    )
    terminal_focus.restore("0x1234")
    assert calls and calls[0][:2] == ["xdotool", "windowactivate"]


def test_linux_capture_swallows_errors_when_xdotool_is_missing(monkeypatch):
    monkeypatch.setattr(terminal_focus.sys, "platform", "linux")

    def raises(*args, **kwargs):
        raise FileNotFoundError("xdotool not found")

    monkeypatch.setattr(terminal_focus.subprocess, "run", raises)
    assert terminal_focus.capture() is None
