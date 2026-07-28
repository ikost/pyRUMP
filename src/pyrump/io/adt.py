r"""Non-Rutherford cross-section tables: ``.adt``, R33 and DSIR 33a.

Port of ``ResRead`` (reswork.c:194-515).

Where the cross-section is not Rutherford -- resonances, light-element EBS --
RUMP reads a measured table and interpolates it **piecewise linearly**. The
manual is untroubled by that: *"The non-smoothness of this description is
completely swamped by convolution with any detector resolution."*

Three dialects, auto-detected from the first non-comment line:

``EARLY_RUMP``
    RUMP's own pre-1999 format. A bare header line ``z1 m1 z2 m2 phi npt``
    followed by ``energy sigma`` pairs. Units are always barns.

``R33``
    The IBANDL/IAEA exchange format. Keyword headers, four data columns
    ``E dE sigma dsigma`` (the uncertainties are read and discarded).

``DSIR 33a``
    An older variant of the same, with two data columns.

Detection is positional, not by extension: a ``VERSION:`` line selects R33 or
DSIR 33a, a leading ``COMMENT:`` starts a comment block, and anything else is
assumed to be an EARLY_RUMP header.

.. note::
   ``phi`` is stored as RUMP's supplement convention. R33 files give
   ``THETA:`` as the true scattering angle, and the reader converts with
   ``phi = 180 - theta`` (reswork.c:361).

.. note::
   Non-zero-Q reactions are rejected outright (reswork.c:355): *"At moment,
   RUMP does not accept non-zero Q reactions"*. So is any exit channel that is
   not the incident particle -- neutrons, gammas and X-rays are refused by
   :func:`parse_reaction`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import numpy as np


class Dialect(Enum):
    EARLY_RUMP = "early_rump"
    R33 = "r33"
    DSIR_33A = "dsir 33a"


class SigmaMode(Enum):
    """How to interpret the tabulated values."""

    BARNS = "barns"
    RELATIVE = "relative"
    """A ratio to Rutherford; multiply by the Rutherford cross-section."""


class AdtError(ValueError):
    """Malformed cross-section table."""


#: ``UNITS:`` keywords RUMP accepts (reswork.c:325-334).
_UNITS = {
    "b/sr": (SigmaMode.BARNS, 1.0),
    "mb/sr": (SigmaMode.BARNS, 0.001),
    "rtr": (SigmaMode.RELATIVE, 1.0),
    "rr": (SigmaMode.RELATIVE, 1.0),
    "relative": (SigmaMode.RELATIVE, 1.0),
}

#: Header keywords that carry no information RUMP uses.
#: MASSES and ZEDS are deliberately among these. RUMP marks both Q_IGNORE
#: (reswork.c:131-132) and takes the nuclide identities from REACTION: alone --
#: the two are redundant and the shipped files disagree on field order
#: (boron.adt gives "11 4 4 11", car_pp.adt "1.0078, 12, 1.0078, 12.0").
_IGNORED = {
    "source", "name", "composition", "zeds", "masses", "subfile", "x4number",
    "serial number", "enddata",
} | {f"address{i}" for i in range(1, 7)}

#: Index is Z, so entry 0 is a placeholder. A neutron is deliberately NOT
#: listed: "n" and nitrogen's "N" differ only in case, and a case-insensitive
#: lookup would silently resolve nitrogen to Z=0. Neutron channels are rejected
#: by _FORBIDDEN before they reach here anyway.
_ELEMENTS = (
    "_ H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te "
    "I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir "
    "Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U"
).split()

#: ``11B(a,a)11B`` -- target(projectile,ejectile)residual.
_REACTION = re.compile(
    r"^\s*(\d*)\s*([A-Za-z]+)\s*\(\s*([^,]+?)\s*,\s*([^)]+?)\s*\)\s*(\d*)\s*([A-Za-z]+)\s*$"
)

#: Shorthands used inside the parentheses (reswork.c:160-176).
_PARTICLES = {
    "a": (2, 4), "alpha": (2, 4),
    "p": (1, 1), "proton": (1, 1),
    "d": (1, 2), "deuteron": (1, 2),
    "t": (1, 3), "triton": (1, 3),
    "h": (2, 3), "he3": (2, 3),
}

#: Exit channels RUMP cannot simulate (reswork.c:172-176).
_FORBIDDEN = {"n", "g", "gamma", "x"}


def _symbol_z(symbol: str) -> int:
    for z, name in enumerate(_ELEMENTS):
        if z and name.lower() == symbol.lower():
            return z
    raise AdtError(f"unknown element symbol {symbol!r}")


def parse_reaction(text: str) -> tuple[int, int, int, int]:
    """Parse ``11B(a,a)11B`` into ``(z1, m1, z2, m2)``.

    Only elastic channels are usable: RUMP requires the ejectile to be the
    projectile and the residual to be the target (reswork.c:186).
    """
    match = _REACTION.match(text)
    if not match:
        raise AdtError(f"cannot parse reaction {text!r}")

    target_mass, target_symbol, incoming, outgoing, residual_mass, residual_symbol = (
        match.groups()
    )

    for particle in (incoming, outgoing):
        if particle.strip().lower() in _FORBIDDEN:
            raise AdtError(
                f"reaction {text!r} has a {particle} channel; RUMP simulates "
                "only charged-particle elastic scattering"
            )

    def resolve(token: str) -> tuple[int, int]:
        key = token.strip().lower()
        # R33 appends the residual's excitation state: p0 is a proton leaving
        # the nucleus in its ground state. Only elastic (state 0) is usable, and
        # the digit carries no mass information.
        state = re.match(r"^([a-z]+)(\d+)$", key)
        if state and state.group(1) in _PARTICLES:
            if state.group(2) != "0":
                raise AdtError(
                    f"{token!r} leaves the nucleus excited; RUMP simulates only "
                    "elastic scattering"
                )
            key = state.group(1)
        if key in _PARTICLES:
            return _PARTICLES[key]
        pair = re.match(r"^(\d+)\s*([A-Za-z]+)$", token.strip())
        if pair:
            return _symbol_z(pair.group(2)), int(pair.group(1))
        raise AdtError(f"unknown particle {token!r}")

    z1, m1 = resolve(incoming)
    z_out, m_out = resolve(outgoing)
    z2 = _symbol_z(target_symbol)
    m2 = int(target_mass) if target_mass else 0

    if (z1, m1) != (z_out, m_out):
        raise AdtError(
            f"reaction {text!r} is not elastic; RUMP requires the ejectile to be "
            "the projectile"
        )
    if _symbol_z(residual_symbol) != z2:
        raise AdtError(f"reaction {text!r} does not return the target nucleus")

    return z1, m1, z2, m2


@dataclass(slots=True)
class CrossSectionTable:
    """One measured excitation function at a fixed angle."""

    energy_keV: np.ndarray
    sigma: np.ndarray
    """Barns/sr, or a ratio to Rutherford when :attr:`mode` is RELATIVE."""

    z1: int = 2
    m1: int = 4
    z2: int = 8
    m2: int = 16
    phi: float = 12.0
    """RUMP's supplement convention: 180 - true scattering angle."""

    mode: SigmaMode = SigmaMode.BARNS
    dialect: Dialect = Dialect.EARLY_RUMP
    comments: list[str] = field(default_factory=list)
    source: Path | None = None

    @property
    def scattering_angle(self) -> float:
        return 180.0 - self.phi

    @property
    def npt(self) -> int:
        return int(self.energy_keV.size)

    def slopes(self) -> np.ndarray:
        """Precomputed piecewise-linear slopes, as RUMP stores (reswork.c:481)."""
        slope = np.zeros_like(self.sigma)
        if self.npt > 1:
            slope[:-1] = np.diff(self.sigma) / np.diff(self.energy_keV)
        return slope

    def __call__(self, energy_keV) -> np.ndarray:
        """Interpolate. Outside the table the endpoints are held flat.

        RUMP instead falls back to Rutherford below the first point and warns on
        overrun; that decision belongs to the caller, which knows the Rutherford
        cross-section. :meth:`covers` reports the valid range.
        """
        energy = np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))
        return np.interp(energy, self.energy_keV, self.sigma)

    def covers(self, energy_keV) -> np.ndarray:
        energy = np.atleast_1d(np.asarray(energy_keV, dtype=np.float64))
        return (energy >= self.energy_keV[0]) & (energy <= self.energy_keV[-1])


def _strip_comment(line: str) -> tuple[str, bool]:
    """Return ``(content, is_comment)``; ``#`` and ``/*`` both introduce one."""
    if line.startswith("#"):
        return line[1:].strip(), True
    if line.startswith("/*"):
        return line[2:].strip(), True
    return line.strip(), False


def read_adt(path: str | Path) -> list[CrossSectionTable]:
    """Read a cross-section file, returning every angle block it contains.

    A single file may hold several blocks separated by ``EndData:``; RUMP
    recycles into header mode and keeps reading (reswork.c:398-402).
    """
    source = Path(path)
    lines = source.read_text(errors="replace").splitlines()

    tables: list[CrossSectionTable] = []
    state = "prescan"
    dialect = Dialect.EARLY_RUMP

    # Defaults are 4He on 16O at 168 degrees (reswork.c:218-228).
    z1, m1, z2, m2 = 2, 4, 8, 16
    phi = 180.0 - 168.0
    mode, unit_scale = SigmaMode.BARNS, 0.001  # mb/sr is the default
    sigma_factor, energy_scale, energy_offset = 1.0, 1.0, 0.0
    comments: list[str] = []
    energies: list[float] = []
    sigmas: list[float] = []
    limit = 65536

    def flush() -> None:
        nonlocal energies, sigmas
        if energies:
            tables.append(
                CrossSectionTable(
                    energy_keV=np.array(energies, dtype=np.float64),
                    sigma=np.array(sigmas, dtype=np.float64),
                    z1=z1, m1=m1, z2=z2, m2=m2, phi=phi,
                    mode=mode, dialect=dialect,
                    comments=list(comments), source=source,
                )
            )
        energies, sigmas = [], []

    for raw in lines:
        content, is_comment = _strip_comment(raw)

        if state == "comment":
            if not content:
                state = "headers"
            elif content:
                comments.append(content)
            continue

        if state == "prescan":
            if content.lower().startswith("version:"):
                value = content.split(":", 1)[1].strip().lower()
                if value.startswith("dsir 33a"):
                    dialect = Dialect.DSIR_33A
                elif value.startswith("r33"):
                    dialect = Dialect.R33
                else:
                    raise AdtError(f"unrecognised version {value!r}")
                state = "headers"
                continue
            if content.lower().startswith("comment:"):
                rest = content.split(":", 1)[1].strip()
                if rest:
                    comments.append(rest)
                state = "comment"
                continue
            if is_comment or not content:
                continue
            # EARLY_RUMP: a bare "z1 m1 z2 m2 phi npt" header.
            fields = content.split()
            if len(fields) < 6:
                raise AdtError(
                    "compatibility-mode header must be 'Z1 M1 Z2 M2 Phi NPT', "
                    f"got {content!r}"
                )
            z1 = int(float(fields[0]))
            m1 = int(float(fields[1]) + 0.5)
            z2 = int(float(fields[2]))
            m2 = int(float(fields[3]) + 0.5)
            phi = float(fields[4])
            limit = int(float(fields[5]))
            mode, unit_scale = SigmaMode.BARNS, 1.0  # always barns
            dialect = Dialect.EARLY_RUMP
            state = "data"
            continue

        if state == "headers":
            if is_comment or not content or ":" not in content:
                continue
            keyword, _, value = content.partition(":")
            keyword, value = keyword.strip().lower(), value.strip()

            if keyword == "comment":
                if value:
                    comments.append(value)
                state = "comment"
            elif keyword == "version":
                low = value.lower()
                dialect = Dialect.DSIR_33A if low.startswith("dsir") else Dialect.R33
            elif keyword == "units":
                if low := _UNITS.get(value.lower()):
                    mode, unit_scale = low
                else:
                    raise AdtError(f"invalid UNITS: {value!r}")
            elif keyword == "reaction":
                z1, m1, z2, m2 = parse_reaction(value)
            elif keyword == "distribution":
                if value.lower() != "energy":
                    raise AdtError("only DISTRIBUTION: ENERGY files are usable")
            elif keyword == "qvalue":
                # The C uses atof(), which consumes only the leading number;
                # several shipped files supply a comma-separated list.
                head = value.replace(",", " ").split()
                if float(head[0]) if head else 0.0:
                    raise AdtError(
                        "non-zero Q reactions are not supported (reswork.c:355)"
                    )
            elif keyword == "theta":
                phi = 180.0 - float(value)
            elif keyword == "enfactors":
                parts = [float(v) for v in value.replace(",", " ").split()]
                energy_scale = parts[0] if parts else 1.0
                energy_offset = parts[1] if len(parts) > 1 else 0.0
            elif keyword == "sigfactors":
                parts = [float(v) for v in value.replace(",", " ").split()]
                sigma_factor = parts[0] if parts else 1.0
            elif keyword == "data":
                limit = 65536
                state = "data"
            elif keyword in _IGNORED:
                pass
            continue

        if state == "data":
            lowered = content.lower()
            if lowered.startswith(("enddata:", "end_data:", "data_end:")):
                flush()
                state = "headers"
                continue
            if is_comment or not content:
                continue

            fields = content.split()
            try:
                if dialect is Dialect.R33:
                    if len(fields) < 4:
                        raise ValueError
                    energy, sigma = float(fields[0]), float(fields[2])
                else:
                    if len(fields) < 2:
                        raise ValueError
                    energy, sigma = float(fields[0]), float(fields[1])
            except ValueError:
                raise AdtError(f"invalid data line: {content!r}") from None

            energy = energy_scale * energy + energy_offset
            sigma = sigma_factor * sigma * unit_scale

            # RUMP drops out-of-order points to keep the table strictly sorted,
            # since the interpolation slopes assume it (reswork.c:427-432).
            if energies and energy <= energies[-1]:
                continue
            energies.append(energy)
            sigmas.append(sigma)
            if len(energies) >= limit:
                flush()
                state = "headers"

    flush()
    if not tables:
        raise AdtError(f"{source}: no cross-section data found")
    return tables


def read_adt_single(path: str | Path) -> CrossSectionTable:
    """Convenience wrapper for the common single-block file."""
    tables = read_adt(path)
    if len(tables) > 1:
        raise AdtError(
            f"{path} holds {len(tables)} angle blocks; use read_adt() instead"
        )
    return tables[0]
