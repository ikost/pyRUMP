"""Give focus back to the terminal after a plot window steals it.

Every GUI backend raises its window when a figure is first shown; on macOS in
particular that makes the new window "key", pulling focus away from whatever
terminal the user was typing in. :func:`capture` records what had focus right
before that happens; :func:`restore` asks the OS to hand it back afterwards.

This is inherently platform-specific and best-effort -- a failure here must
never break a plotting command, so every backend swallows its own errors.
"""

from __future__ import annotations

import subprocess
import sys

_TIMEOUT = 1.0


def capture() -> object | None:
    """A token identifying whatever has focus right now, or ``None``."""
    if sys.platform == "darwin":
        return _capture_macos()
    if sys.platform == "win32":
        return _capture_windows()
    if sys.platform.startswith("linux"):
        return _capture_linux()
    return None


def restore(token: object | None) -> None:
    """Best-effort: give focus back to whatever :func:`capture` recorded."""
    if token is None:
        return
    if sys.platform == "darwin":
        _restore_macos(token)
    elif sys.platform == "win32":
        _restore_windows(token)
    elif sys.platform.startswith("linux"):
        _restore_linux(token)


# -- macOS: osascript, no extra dependency -----------------------------------
#
# Re-activating another app requires driving "System Events", which macOS
# gates behind a one-time Automation permission prompt the first time this
# runs -- there is no way around that short of a compiled helper.


def _capture_macos() -> str | None:
    try:
        result = subprocess.run(
            [
                "osascript", "-e",
                'tell application "System Events" to get name of first '
                "process whose frontmost is true",
            ],
            capture_output=True, text=True, timeout=_TIMEOUT,
        )
        name = result.stdout.strip()
        return name or None
    except Exception:
        return None


def _restore_macos(app_name: str) -> None:
    try:
        subprocess.run(
            ["osascript", "-e", f'tell application "{app_name}" to activate'],
            capture_output=True, timeout=_TIMEOUT,
        )
    except Exception:
        pass


# -- Windows: ctypes + user32, no extra dependency ---------------------------


def _capture_windows() -> int | None:
    try:
        import ctypes

        hwnd = ctypes.windll.user32.GetForegroundWindow()
        return hwnd or None
    except Exception:
        return None


def _restore_windows(hwnd: int) -> None:
    try:
        import ctypes

        ctypes.windll.user32.SetForegroundWindow(hwnd)
    except Exception:
        pass


# -- Linux: xdotool if it happens to be installed, else a no-op -------------
#
# X11-only, and only as good as the window manager's cooperation; there is no
# analogous trick on Wayland, so this quietly does nothing there.


def _capture_linux() -> str | None:
    try:
        result = subprocess.run(
            ["xdotool", "getactivewindow"],
            capture_output=True, text=True, timeout=_TIMEOUT,
        )
        window_id = result.stdout.strip()
        return window_id or None
    except Exception:
        return None


def _restore_linux(window_id: str) -> None:
    try:
        subprocess.run(
            ["xdotool", "windowactivate", window_id],
            capture_output=True, timeout=_TIMEOUT,
        )
    except Exception:
        pass
