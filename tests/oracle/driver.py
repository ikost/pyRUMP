"""Drive the legacy RUMP C binary as a numerical reference oracle.

The shipped ``C-code/bin/rump`` is a curses-style TTY application: it produces no
output at all when stdin is a pipe, so it must be driven through a pseudo-terminal.
This module wraps that in a batch interface.

The C tree is not redistributed with pyRUMP. Set ``PYRUMP_C_REFERENCE`` to its root
(the directory containing ``bin/rump``), or place it at ``C-code/`` beside the repo.
Every entry point degrades to ``None``/skip when the tree is absent.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

# RUMP writes VT100 escapes and cursor addressing throughout; strip for parsing.
_ANSI = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]|\x1b[()][A-Z0-9]|\r")

# Prompts RUMP uses at the top level. It rotates through these jokingly
# ("Your wish?", "Yes Master?", ...), so we match the trailing "? " generically
# as well as the SIM/PERT sub-shell prompts.
_PROMPT = re.compile(
    r"(?:SIM Command|PERT Command|CONFIG Command):\s*$"
    r"|(?:Your wish|Yes Master|Yes dear|Here I am|Ready if you are|Up periscope)[?!]\s*$"
)


def reference_root() -> Path | None:
    """Locate the legacy C tree, or return None if unavailable."""
    env = os.environ.get("PYRUMP_C_REFERENCE")
    candidates = [Path(env)] if env else []
    candidates.append(Path(__file__).resolve().parents[2] / "C-code")
    for root in candidates:
        if (root / "bin" / "rump").is_file():
            return root
    return None


def available() -> bool:
    return reference_root() is not None


@dataclass
class OracleResult:
    """Outcome of one batch RUMP session."""

    transcript: str
    workdir: Path

    def output(self, name: str) -> Path:
        return self.workdir / name


class RumpOracle:
    """Batch driver for the legacy RUMP binary.

    RUMP resolves its data tables (atom4.dat, pscoef.dat, newstop.kal, density.tab)
    relative to the binary's directory, so the child runs with cwd=bin/. Output files
    are therefore written with absolute paths into a scratch directory.
    """

    def __init__(self, root: Path | None = None, timeout: float = 60.0):
        resolved = root or reference_root()
        if resolved is None:
            raise RuntimeError(
                "Legacy RUMP C tree not found. Set PYRUMP_C_REFERENCE to the directory "
                "containing bin/rump."
            )
        self.root = resolved
        self.binary = resolved / "bin" / "rump"
        self.timeout = timeout

    def run(self, commands: list[str], workdir: Path | None = None) -> OracleResult:
        """Execute `commands` in a fresh RUMP session and return the transcript."""
        import pexpect  # imported lazily so the module imports without the dev extra

        work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="rump-oracle-"))
        work.mkdir(parents=True, exist_ok=True)

        child = pexpect.spawn(
            str(self.binary),
            cwd=str(self.binary.parent),
            timeout=self.timeout,
            encoding="utf-8",
            codec_errors="replace",
            env={**os.environ, "TERM": "vt100"},
            dimensions=(24, 80),
        )
        transcript: list[str] = []
        try:
            # Startup: splash screen (with a deliberate pause), data-table loading,
            # then the rump.ini macro. Needs a much longer quiet window than a
            # normal command before we can conclude RUMP is idle at a prompt.
            self._settle(child, transcript, quiet_for=3.0)
            for cmd in commands:
                child.send(cmd + "\n")
                self._settle(child, transcript)
        finally:
            try:
                child.send("quit\n")
                child.expect([pexpect.EOF, pexpect.TIMEOUT], timeout=5)
            except Exception:
                pass
            child.close(force=True)

        return OracleResult(transcript=_ANSI.sub("", "".join(transcript)), workdir=work)

    def _settle(self, child, transcript: list[str], quiet_for: float = 0.4) -> None:
        """Read until RUMP stops emitting, i.e. it is waiting at a prompt."""
        import pexpect

        while True:
            try:
                child.expect(r".+", timeout=quiet_for, searchwindowsize=None)
                transcript.append(child.after)
            except pexpect.TIMEOUT:
                return
            except pexpect.EOF:
                if child.before:
                    transcript.append(child.before)
                return

    def simulate_to_ascii(
        self,
        macro: str,
        out: Path,
        *,
        data_file: Path | str,
        setup: list[str] | None = None,
    ) -> Path:
        """Run a SIM macro against a loaded data buffer and write the theory spectrum.

        `data_file` is required: RUMP simulates onto the *active* buffer's channel
        grid and beam parameters (`RbsCopySpectrum(ALTBUF, ibuf)`, creatr.c:283), so
        with no data loaded the simulation produces npt=0. Use one of the shipped
        `data/Fixed/*.rbs` fixtures, or any real spectrum.

        `macro` is the body of SIM commands, without the enclosing `sim`/`return`.
        """
        out = Path(out).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.unlink(missing_ok=True)

        cmds = [f"read '{Path(data_file).resolve()}'"]
        cmds.extend(setup or [])
        cmds.append("sim")
        cmds.extend(line.strip() for line in macro.strip().splitlines() if line.strip())
        # NOTE: the SIM keyword is misspelled "recalculculate" in the C command table
        # (sim2.c:252) with a 5-character minimum, so "recalculate" does NOT match and
        # silently falls through. "recal" is the only safe spelling.
        cmds.append("recal")
        cmds.append("return")
        cmds.append("pointat theory")
        cmds.append(f"wrascii '{out}'")

        result = self.run(cmds, workdir=out.parent)
        if not out.exists():
            raise RuntimeError(f"RUMP produced no output.\nTranscript:\n{result.transcript[-2000:]}")
        return out


def read_wrascii(path: Path) -> tuple[dict[str, str], list[float]]:
    """Parse RUMP's ``wrascii`` output into (metadata, counts).

    Layout is a keyword header block, the literal line ``Swallow``, then one count
    per line::

        Spectrum    RBS
        Ident      'Simulation of Si'
        Charge      10.000000     MeV  3.000000
        Conversion  7.815000 65.646000
        ...
        Swallow
        0.000000
        ...
    """
    lines = Path(path).read_text().splitlines()
    try:
        split = next(i for i, ln in enumerate(lines) if ln.strip().lower() == "swallow")
    except StopIteration as exc:
        raise ValueError(f"{path}: no 'Swallow' marker; not a wrascii file") from exc

    meta: dict[str, str] = {}
    for ln in lines[:split]:
        parts = ln.split(None, 1)
        if len(parts) == 2:
            meta[parts[0]] = parts[1].strip()

    counts: list[float] = []
    for ln in lines[split + 1 :]:
        ln = ln.strip()
        if ln:
            counts.append(float(ln))
    return meta, counts


__all__ = [
    "RumpOracle",
    "OracleResult",
    "available",
    "reference_root",
    "read_wrascii",
]
