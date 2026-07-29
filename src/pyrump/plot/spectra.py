"""matplotlib helpers for RBS spectra.

pyRUMP does not port genplot — these are ordinary matplotlib figures, returned
so callers can restyle or embed them.

Conventions chosen for readability rather than to imitate RUMP's output:
measured data as a step (it is binned counts, not a continuous curve),
simulation as a line over it, and residuals in units of sigma against a shaded
+/-1 band.
"""

from __future__ import annotations

import numpy as np

from ..fit.objective import poisson_residuals
from ..model.spectrum import Spectrum


def _axes(ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))
    return ax


def plot_spectrum(
    spectrum: Spectrum,
    *,
    ax=None,
    energy_axis: bool = True,
    label: str | None = None,
    **kwargs,
):
    """Draw one spectrum against energy (default) or channel number."""
    ax = _axes(ax)
    x = spectrum.energies if energy_axis else np.arange(spectrum.counts.size)
    ax.plot(x, spectrum.counts, label=label, **kwargs)
    ax.set_xlabel("Energy (keV)" if energy_axis else "Channel")
    ax.set_ylabel("Counts")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(bottom=0)
    if label:
        ax.legend(frameon=False)
    return ax


def plot_comparison(
    data: Spectrum,
    simulation: Spectrum,
    *,
    energy_axis: bool = True,
    residuals: bool = True,
    window: np.ndarray | None = None,
    figsize: tuple[float, float] = (9, 6),
):
    """Measured data with a simulation over it, and optionally residuals.

    Residuals are the Poisson residuals the fit actually minimises, not
    ``data - model``, so what is plotted is what was optimised. Channels where
    the model predicts zero are omitted: they contribute nothing and plotting
    them as zero would imply agreement.
    """
    import matplotlib.pyplot as plt

    n = min(data.counts.size, simulation.counts.size)
    observed, expected = data.counts[:n], simulation.counts[:n]
    x = data.energies[:n] if energy_axis else np.arange(n)

    if residuals:
        figure, (top, bottom) = plt.subplots(
            2, 1, figsize=figsize, sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        )
    else:
        figure, top = plt.subplots(figsize=figsize)
        bottom = None

    top.step(x, observed, where="mid", lw=0.9, color="0.35", label="data")
    top.plot(x, expected, lw=1.6, color="crimson", label="simulation")
    top.set_ylabel("Counts")
    top.set_ylim(bottom=0)
    top.legend(frameon=False)

    if window is not None:
        mask = np.asarray(window, dtype=bool)[:n]
        for start, stop in _runs(mask):
            top.axvspan(x[start], x[stop], color="steelblue", alpha=0.10, lw=0)

    if bottom is not None:
        values, _ = poisson_residuals(observed, expected)
        shown = expected > 0
        bottom.axhspan(-1, 1, color="0.85", zorder=0)
        bottom.axhline(0, color="0.4", lw=0.8)
        bottom.plot(x[shown], values[shown], lw=0.9, color="crimson")
        bottom.set_ylabel(r"residual ($\sigma$)")
        bottom.set_xlabel("Energy (keV)" if energy_axis else "Channel")
        limit = max(3.0, float(np.abs(values[shown]).max()) if shown.any() else 3.0)
        bottom.set_ylim(-limit, limit)
    else:
        top.set_xlabel("Energy (keV)" if energy_axis else "Channel")

    top.set_xlim(x.min(), x.max())
    return figure


def plot_depth_profile(
    grid, element_symbols: list[str], *, ax=None, atomic_fraction: bool = True
):
    """Composition against depth, one line per element."""
    ax = _axes(ax)
    depth = np.concatenate([[0.0], grid.depth])
    for column, symbol in enumerate(element_symbols):
        values = grid.composition[:, column]
        if atomic_fraction:
            with np.errstate(divide="ignore", invalid="ignore"):
                values = np.where(grid.areal_density > 0, values / grid.areal_density, 0)
        # Step through the slabs: composition is constant within each.
        ax.step(depth, np.concatenate([[values[0]], values]), where="pre", label=symbol)

    ax.set_xlabel(r"Depth (10$^{15}$ atoms/cm$^2$)")
    ax.set_ylabel("Atomic fraction" if atomic_fraction else r"10$^{15}$ atoms/cm$^2$")
    ax.set_xlim(0, depth[-1])
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    return ax


def _runs(mask: np.ndarray):
    """Yield (start, stop) index pairs for each contiguous True run."""
    if not mask.any():
        return
    edges = np.diff(mask.astype(int))
    starts = list(np.flatnonzero(edges == 1) + 1)
    stops = list(np.flatnonzero(edges == -1))
    if mask[0]:
        starts.insert(0, 0)
    if mask[-1]:
        stops.append(len(mask) - 1)
    yield from zip(starts, stops)
