"""Buffers, plot state, and the session that owns them.

RUMP keeps spectra in numbered *buffers* (concepts.htm): ``-1`` scratch, ``0``
the simulation, ``1..10`` experimental data, with exactly one designated ACTIVE
and operated on implicitly by most commands. pyRUMP's library has no such
concept -- a :class:`~pyrump.model.spectrum.Spectrum` is only counts plus a
calibration -- so the buffer is the one genuinely new abstraction the shell adds.

A :class:`Buffer` is the spectrum together with the experimental parameters
needed to analyse it, which is very nearly
:class:`~pyrump.fit.parameters.FitInputs` and exactly what
:class:`~pyrump.io.rbs.RbsSpectrum` already carries.

Two deliberate divergences from the original:

* **No fixed ring.** RUMP recycled ten slots destructively ("the original buffer
  10 is lost in the process"). That was a memory constraint, not a feature.
  Buffers here are a list that grows, still addressed by number.
* **Slot 0 is still the simulation**, and keeps RUMP's best idea: there is no
  "simulate" command. Buffer 0 is recomputed lazily whenever the sample
  description or the active buffer's parameters change.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from ..model.detector import Measurement
from ..model.geometry import Geometry
from ..model.spectrum import Calibration, Spectrum
from ..script.lcm import Script
from ..sim.engine import Beam

#: Names RUMP accepts in place of a buffer number (cmds.htm, PLOT).
SPECIAL_NAMES = {
    "sim": 0,
    "alt": 0,
    "theory": 0,
    "main": 1,
}


@dataclass(slots=True)
class Buffer:
    """One spectrum with the parameters needed to interpret it."""

    spectrum: Spectrum
    beam: Beam = field(default_factory=Beam)
    geometry: Geometry = field(default_factory=lambda: Geometry(theta=0.0, phi=10.0))
    measurement: Measurement = field(default_factory=Measurement)

    name: str = ""
    path: Path | None = None
    identifier: str = ""
    date: str = ""
    livetime: str = ""
    comments: list[str] = field(default_factory=list)

    @property
    def calibration(self) -> Calibration:
        return self.spectrum.calibration

    @property
    def n_channels(self) -> int:
        return int(self.spectrum.counts.size)

    @classmethod
    def from_rbs(cls, source, path: Path | None = None) -> "Buffer":
        """Build a buffer from :func:`pyrump.io.rbs.read_rbs`."""
        return cls(
            spectrum=source.to_spectrum(),
            beam=Beam(
                e0_MeV=source.e0_MeV,
                z=source.zbeam or 2,
                mass=source.mbeam or 4.0026,
                charge_state=source.measurement.charge_state or 1,
            ),
            geometry=source.geometry,
            measurement=source.measurement,
            name=path.name if path else source.identifier,
            path=path,
            identifier=source.identifier,
            date=source.date,
            livetime=source.livetime,
            comments=list(source.comments),
        )

    @classmethod
    def from_ascii(cls, source, path: Path | None = None) -> "Buffer":
        """Build a buffer from :func:`pyrump.io.ascii.read_ascii`.

        ASCII files carry no metadata, so everything but the counts takes its
        default and the user is expected to set it with the metadata commands.
        """
        counts = np.asarray(source.counts, dtype=float)
        calibration = getattr(source, "calibration", None) or Calibration(
            npt=counts.size
        )
        return cls(
            spectrum=Spectrum(counts=counts, calibration=calibration),
            name=path.name if path else "",
            path=path,
            identifier=getattr(source, "identifier", ""),
        )

    def to_rbs(self):
        """Convert back for :func:`pyrump.io.rbs.write_rbs`."""
        from ..io.rbs import RbsSpectrum

        return RbsSpectrum(
            counts=self.spectrum.counts,
            calibration=self.calibration,
            geometry=self.geometry,
            measurement=self.measurement,
            e0_MeV=self.beam.e0_MeV,
            zbeam=self.beam.z,
            mbeam=self.beam.mass,
            identifier=self.identifier,
            date=self.date,
            livetime=self.livetime,
            comments=list(self.comments),
        )

    def copy(self) -> "Buffer":
        """A deep-enough copy: the frozen parameter objects can be shared."""
        return replace(
            self,
            spectrum=Spectrum(
                counts=self.spectrum.counts.copy(), calibration=self.calibration
            ),
            comments=list(self.comments),
        )

    # -- metadata setters -------------------------------------------------
    # Calibration, Geometry and Measurement are frozen, so every setter
    # rebuilds rather than mutates -- the same discipline fit/parameters.py uses.

    def set_calibration(self, **changes) -> None:
        self.spectrum.calibration = replace(self.calibration, **changes)

    def set_geometry(self, **changes) -> None:
        self.geometry = replace(self.geometry, **changes)

    def set_measurement(self, **changes) -> None:
        self.measurement = replace(self.measurement, **changes)

    def display_path(self) -> str:
        """Path for listings: relative to the working directory when under it.

        Buffers store absolute paths so that lookups survive a ``CD``; this
        keeps the printed form short anyway.
        """
        if self.path is None:
            return "(none)"
        try:
            return str(self.path.relative_to(Path.cwd()))
        except ValueError:
            return str(self.path)

    def describe(self) -> str:
        """The ``ACTIVE`` listing: everything RUMP shows for one buffer."""
        c, g, m, b = self.calibration, self.geometry, self.measurement, self.beam
        lines = [
            f"  File       {self.display_path()}",
            f"  Identifier {self.identifier}",
            f"  Date       {self.date}    Livetime {self.livetime}",
            f"  Beam       Z={b.z} mass={b.mass:.4f}  {b.e0_MeV:.4f} MeV"
            f"  charge state {m.charge_state}",
            f"  Geometry   {g.kind.name.lower()}  theta {g.theta:g}  phi {g.phi:g}"
            f"  psi {g.psi:g}   (scattering angle {g.scattering_angle:g})",
            f"  Conversion {c.kevch:g} keV/ch   offset {c.kev0:g} keV"
            f"   first channel {c.first:g}   {self.n_channels} channels",
            f"  Detector   FWHM {m.fwhm_keV:g} keV   omega {m.omega_msr:g} msr"
            f"   tau {m.tau_us:g} us",
            f"  Dose       {m.charge_uC:g} uC   corr {m.correction:g}"
            f"   current {m.current_nA:g} nA",
            f"  Total      {self.spectrum.total():.1f} counts",
        ]
        return "\n".join(lines)


@dataclass(slots=True)
class BufferSet:
    """The numbered buffers, and which one is ACTIVE.

    Index 0 is always the simulation. Data buffers start at 1. Unlike the
    original there is no upper limit and nothing is silently destroyed.
    """

    slots: list[Buffer | None] = field(default_factory=lambda: [None])
    active: int = 0

    def __len__(self) -> int:
        return len(self.slots)

    def __getitem__(self, index: int) -> Buffer:
        buffer = self.get(index)
        if buffer is None:
            raise KeyError(f"buffer {index} is empty")
        return buffer

    def get(self, index: int) -> Buffer | None:
        if 0 <= index < len(self.slots):
            return self.slots[index]
        return None

    def set(self, index: int, buffer: Buffer | None) -> None:
        while len(self.slots) <= index:
            self.slots.append(None)
        self.slots[index] = buffer

    @property
    def active_buffer(self) -> Buffer | None:
        return self.get(self.active)

    def require_active(self) -> Buffer:
        buffer = self.active_buffer
        if buffer is None:
            raise KeyError("no active buffer: read a spectrum first (GET <file>)")
        return buffer

    def first_free(self) -> int:
        """Lowest empty data slot, appending one if all are full."""
        for index in range(1, len(self.slots)):
            if self.slots[index] is None:
                return index
        self.slots.append(None)
        return len(self.slots) - 1

    def find_path(self, path: Path) -> int | None:
        """The buffer already holding ``path``, if any.

        RUMP checks the full pathname before reading a file again (cmds.htm,
        PLOT), so re-plotting a file just re-selects its buffer.

        Both sides are absolute: buffers store a resolved path at load time, so
        that a ``CD`` between two ``GET``s of the same file cannot make it look
        like two different files and load it twice.
        """
        resolved = path.resolve()
        for index, buffer in enumerate(self.slots):
            if buffer is not None and buffer.path is not None:
                if buffer.path == resolved:
                    return index
        return None

    def load(self, buffer: Buffer, index: int | None = None) -> int:
        """Place a buffer, defaulting to the lowest free data slot."""
        if index is None:
            index = self.first_free()
        self.set(index, buffer)
        return index

    def release(self, index: int) -> None:
        if index == 0:
            raise KeyError("buffer 0 is the simulation and cannot be released")
        self.set(index, None)
        if self.active == index:
            self.active = next(
                (i for i, b in enumerate(self.slots) if i and b is not None), 0
            )

    def clear(self) -> None:
        self.slots = [self.slots[0] if self.slots else None]
        self.active = 0

    def listing(self) -> str:
        """The ``BUFFERS`` table."""
        lines = ["  #  ACT  channels  total counts  file"]
        for index, buffer in enumerate(self.slots):
            if buffer is None and index:
                continue
            mark = "*" if index == self.active else " "
            if buffer is None:
                lines.append(f"{index:3d}   {mark}   (empty simulation)")
                continue
            label = buffer.display_path() if buffer.path else (buffer.identifier or "-")
            lines.append(
                f"{index:3d}   {mark}  {buffer.n_channels:8d}"
                f"  {buffer.spectrum.total():12.1f}  {label}"
            )
        return "\n".join(lines)


@dataclass(slots=True)
class XeqFrame:
    """One nested ``XEQ`` file's raw lines, and how far dispatch has gotten.

    :func:`~pyrump.shell.repl.execute_file` advances ``index`` one line at a
    time. ``SWALLOW`` (rump.py) is the one command that needs more than its own
    line -- it reaches into the innermost frame and consumes the following raw
    data lines itself, exactly as RUMP's own macro lexer reads channel data
    straight off the active input stream (bmanip.c: B1_SWALLOW).
    """

    lines: list[str]
    index: int = 0


@dataclass(slots=True)
class PlotState:
    """Everything the plotting commands manipulate.

    Mirrors RUMP's persistent plot parameters (cmds.htm, "Plot Scaling and Axes
    Layout"): the channel region, the yield range, the yield scaling, and
    whether yield is normalised or raw.
    """

    low: int | None = None
    high: int | None = None
    ylow: float | None = None
    yhigh: float | None = None
    yscale: str = "linear"          # linear | sqrt | log
    normalized: bool = False
    energy_axis: bool = False       # RUMP plots against channel by default
    labels: bool = True
    forcex: bool = False
    autoids: bool = False

    def region(self, n_channels: int) -> tuple[int, int]:
        """The channel range to draw, defaulting to the whole spectrum."""
        low = 0 if self.low is None else max(0, self.low)
        high = n_channels - 1 if self.high is None else min(n_channels - 1, self.high)
        if high <= low:
            raise ValueError(f"empty plot region: {low} to {high}")
        return low, high

    def describe(self) -> str:
        """The ``PARMS`` listing."""
        span = (
            f"{self.low if self.low is not None else 'auto'}"
            f" to {self.high if self.high is not None else 'auto'}"
        )
        counts = (
            f"{self.ylow if self.ylow is not None else 'auto'}"
            f" to {self.yhigh if self.yhigh is not None else 'auto'}"
        )
        return "\n".join(
            [
                f"  Region     {span} (channels)",
                f"  Counts     {counts}",
                f"  Yield axis {self.yscale}"
                f", {'normalized' if self.normalized else 'raw'}",
                f"  X axis     {'energy' if self.energy_axis else 'channel'}"
                f"{', forced' if self.forcex else ''}",
                f"  Labels     {'on' if self.labels else 'off'}",
            ]
        )


@dataclass(slots=True)
class Session:
    """Everything the interactive shell owns.

    The atomic/stopping tables are built once here and reused for every
    simulation, as ``docs/usage.md`` requires.
    """

    table: object
    registry: object
    densities: object
    data: Path

    buffers: BufferSet = field(default_factory=BufferSet)
    plot: PlotState = field(default_factory=PlotState)

    #: The SIM sample description.
    script: Script = field(default_factory=Script)
    #: The SIM layer editor over :attr:`script`, holding the layer pointer.
    editor: object | None = None
    #: Set when the sample or the active buffer's parameters change; buffer 0
    #: is recomputed on next use. RUMP has no explicit "simulate" command.
    dirty: bool = True

    #: PERT state, populated by :mod:`pyrump.shell.commands.pert`.
    pert: object | None = None

    figure: object | None = None
    traces: list = field(default_factory=list)
    log_file: object | None = None
    echo: bool = False
    #: XEQ nesting depth, guarded against a macro that calls itself.
    xeq_depth: int = 0
    #: One :class:`XeqFrame` per nested XEQ, innermost last.
    xeq_stack: list[XeqFrame] = field(default_factory=list)
    #: PUSHDIR/POPDIR stack of directories to come back to.
    directory_stack: list[Path] = field(default_factory=list)

    #: INTSET's two mode flags, shared by INTEGRAL/THICKNESS (RbsThickn's
    #: persistent `interp`/`qmode` statics, anlytc.c:1506).
    integration_interp: bool = False
    integration_qmode: int = 0

    def xeq_frame(self) -> XeqFrame | None:
        """The innermost running macro's line frame, or ``None`` outside XEQ."""
        return self.xeq_stack[-1] if self.xeq_stack else None

    @classmethod
    def create(cls, data: str | None = None) -> "Session":
        from ..cli._common import data_dir, load_tables

        directory = data_dir(data)
        table, registry, densities = load_tables(directory)
        return cls(table=table, registry=registry, densities=densities, data=directory)

    # -- buffer resolution -------------------------------------------------

    def resolve(self, token: str) -> int:
        """Turn a buffer token into an index.

        Accepts a number, one of RUMP's special names (``SIM``/``ALT``/
        ``THEORY``/``MAIN``/``NOW``/``LAST``), or a filename already held in a
        buffer. Does **not** read files -- callers that should do so check for
        :class:`KeyError` and read.
        """
        text = token.strip()
        if not text:
            raise KeyError("expected a buffer")
        try:
            return int(text)
        except ValueError:
            pass
        lowered = text.lower()
        if lowered in ("now", "last"):
            return self.buffers.active
        if lowered in SPECIAL_NAMES:
            return SPECIAL_NAMES[lowered]
        found = self.buffers.find_path(Path(text))
        if found is not None:
            return found
        raise KeyError(f"no buffer named {token!r}")

    def touch(self) -> None:
        """Mark the simulation stale."""
        self.dirty = True

    # -- the implicit simulation ------------------------------------------

    def simulation(self) -> Buffer:
        """Buffer 0, recomputed if the sample or active parameters changed.

        This is RUMP's design: ``PLOT 0``/``COMPARE`` update the theory spectrum
        when necessary rather than making the user ask (sim.htm, "computational
        routines within RUMP update the theoretical spectrum as required").
        """
        from ..script.lcm import to_sample
        from ..sim.engine import simulate

        existing = self.buffers.get(0)
        if existing is not None and not self.dirty:
            return existing

        if not self.script.layers:
            raise KeyError("no sample described: use SIM to build one")
        reference = self.buffers.require_active()

        sample = to_sample(self.script, self.table, self.densities)
        spectrum = simulate(
            sample,
            reference.beam,
            reference.geometry,
            self.registry,
            self.table,
            reference.calibration,
            reference.measurement,
        )
        buffer = Buffer(
            spectrum=spectrum,
            beam=reference.beam,
            geometry=reference.geometry,
            measurement=reference.measurement,
            name="SIM",
            identifier=self.script.description or "pyRUMP simulation",
        )
        self.buffers.set(0, buffer)
        self.dirty = False
        return buffer

    def write_log(self, line: str) -> None:
        if self.log_file is not None:
            self.log_file.write(line.rstrip() + "\n")
            self.log_file.flush()
