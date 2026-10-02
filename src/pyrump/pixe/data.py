"""A measured PIXE spectrum, as a buffer carries it next to its RBS one."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from ..model.spectrum import Calibration, Spectrum


@dataclass(slots=True)
class PixeData:
    """One PIXE spectrum with what is needed to interpret it.

    The beam, tilt and charge are not here: they belong to the buffer, which
    shares them with its RBS spectrum (one run, one charge integrator).
    """

    spectrum: Spectrum
    path: Path | None = None
    identifier: str = ""
    date: str = ""
    live_time_s: float | None = None
    real_time_s: float | None = None

    @property
    def calibration(self) -> Calibration:
        return self.spectrum.calibration

    @property
    def live_fraction(self) -> float:
        """LT/RT: the fraction of the run the detector was counting."""
        if self.live_time_s and self.real_time_s:
            return self.live_time_s / self.real_time_s
        return 1.0

    def copy(self) -> "PixeData":
        return replace(
            self,
            spectrum=Spectrum(
                counts=self.spectrum.counts.copy(), calibration=self.calibration
            ),
        )

    def describe(self) -> str:
        """One line for a buffer's ``ACTIVE`` listing."""
        name = self.path.name if self.path else (self.identifier or "-")
        c = self.calibration
        times = (
            f"   live {self.live_time_s:g} s / real {self.real_time_s:g} s"
            if self.live_time_s and self.real_time_s
            else ""
        )
        return (
            f"  PIXE       {name}   {c.npt} channels, {c.kevch * 1000:.4g} eV/ch,"
            f" offset {c.kev0 * 1000:.4g} eV{times}"
        )

    @classmethod
    def read(cls, path: Path, default_calibration: Calibration) -> tuple["PixeData", list[str]]:
        """Read a PIXE spectrum file. Returns the data and notices for the user.

        Only Oxford Instruments ASCII (RC43's ``.PIX``) is known so far.
        Without a usable header calibration, ``default_calibration``'s gain
        and offset are used, on the file's own channel numbers.
        """
        from ..io.oxford import read_oxford

        source = read_oxford(path)
        calibration = source.calibration or replace(
            default_calibration,
            first=float(source.first_channel),
            npt=source.counts.size,
        )
        data = cls(
            spectrum=Spectrum(counts=source.counts, calibration=calibration),
            path=Path(path),
            identifier=source.identifier,
            date=source.date,
            live_time_s=source.live_time_s,
            real_time_s=source.real_time_s,
        )
        return data, list(source.notices)
