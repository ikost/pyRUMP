"""Draw the figures of the manual's Experimental geometry page.

    PYTHONPATH=src python3.14 tools/geometry_figures.py

writes docs/assets/geometry-ibm.png: GEOMETRY IBM seen along the tilt axis --
beam, detector and normal all lie in the scattering plane, so a flat drawing
shows every angle true. The exit angle in each panel's box is computed by
pyRUMP itself (``pyrump.model.geometry.Geometry``), so the figure can't
drift from the code.
"""

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Arc, Polygon  # noqa: E402

THETA, PHI = 25.0, 35.0  # as in the other sketches


def polar(v) -> float:
    return math.degrees(math.atan2(v[1], v[0]))


def vec(deg_from_up: float) -> np.ndarray:
    """A unit vector at ``deg_from_up`` from straight up (towards the beam's
    source), positive to the right."""
    r = math.radians(deg_from_up)
    return np.array([math.sin(r), math.cos(r)])


def arrow(ax, start, end, color, lw=2.5, ls="-"):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": lw, "ls": ls,
                            "mutation_scale": 18})


def angle_arc(ax, a, b, radius, color, label, at, ha="center"):
    """An arc between ``a`` and ``b``, labelled at ``at`` = (degrees from
    up, distance from the origin)."""
    t1, t2 = sorted((polar(a), polar(b)))
    ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=t1, theta2=t2, color=color, lw=2))
    x, y = at[1] * vec(at[0])
    ax.text(x, y, label, color=color, fontsize=11, fontweight="bold", ha=ha, va="center")


def rump_exit(theta):
    from pyrump.model.geometry import Geometry, GeometryKind

    geometry = Geometry(theta=theta, phi=PHI, kind=GeometryKind.IBM)
    return math.degrees(math.acos(1.0 / geometry.sec_out))


def panel(ax, theta, title):
    s = vec(0)               # back towards the beam's source
    d = vec(-PHI)            # the RBS detector, to the left
    n = vec(theta)           # IBM: THETA > 0 turns the normal away from the detector
    psi = math.degrees(math.acos(np.dot(n, d)))

    # The sample, edge-on: a slab perpendicular to the normal, material below.
    along = np.array([n[1], -n[0]])
    corners = [1.0 * along, -1.0 * along, -1.0 * along - 0.12 * n, 1.0 * along - 0.12 * n]
    ax.add_patch(Polygon(corners, closed=True, facecolor="0.75", edgecolor="0.35"))
    ax.text(*(0.95 * along - 0.3 * n), "sample", color="0.3", fontsize=10, ha="center")

    arrow(ax, (0, 2.25), (0, 0.04), "black")
    ax.text(0.06, 2.15, "beam", fontsize=11)
    arrow(ax, (0, 0), 1.95 * d, "darkmagenta")
    ax.text(*(2.02 * d), "RBS detector", color="darkmagenta", fontsize=11, ha="right")
    arrow(ax, (0, 0), 1.45 * n, "royalblue")
    if theta > 0:
        ax.text(*(1.5 * n + [0.0, 0.08]), "normal", color="royalblue", fontsize=11, ha="left")
    else:
        ax.text(*(1.45 * n + [0.08, -0.12]), "normal", color="royalblue", fontsize=11,
                ha="left")
    ax.plot([0, 0], [0, 1.35], color="royalblue", ls="--", lw=1)
    if theta > 0:
        ax.text(-0.05, 1.5, "untilted normal\n(along the beam)", color="royalblue",
                fontsize=8, ha="right", va="top")
    else:
        ax.text(0.04, 1.3, "untilted normal\n(along the beam)", color="royalblue",
                fontsize=8, ha="left", va="top")

    # The tilt axis points out of the page.
    ax.plot(0, 0, marker="o", ms=11, mfc="white", mec="purple", mew=2)
    ax.plot(0, 0, marker=".", ms=6, color="purple")
    ax.text(0.1, -0.32, "tilt axis\n(out of the page)", color="purple", fontsize=9)

    angle_arc(ax, n, s, 0.5, "crimson", f"THETA {abs(theta):g}°", at=(14, 0.62), ha="left")
    angle_arc(ax, s, d, 1.62, "darkorange", f"PHI {PHI:g}°", at=(-17.5, 1.76))
    if theta > 0:
        angle_arc(ax, n, d, 0.85, "darkmagenta", f"PSI {psi:.0f}°", at=(-8, 0.98), ha="right")
    else:
        angle_arc(ax, n, d, 0.75, "darkmagenta", f"PSI {psi:.0f}°", at=(-44, 0.78), ha="right")

    out = rump_exit(theta)
    lines = [
        "RUMP:",
        "  GEOMETRY IBM",
        f"  THETA {theta:g}     beam to the normal",
        f"  PHI {PHI:g}        180° - scattering angle ({180 - PHI:g}°)",
        "  -> detector to the normal:",
        f"     PSI = |THETA + PHI| = |{theta:g} + {PHI:g}| = {out:g}°",
    ]
    ax.text(-2.3, -1.15, "\n".join(lines), family="monospace", fontsize=9.5, va="top",
            bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "0.6"})

    ax.set_title(title, fontsize=12)
    ax.set_xlim(-2.4, 2.0)
    ax.set_ylim(-2.0, 2.4)
    ax.set_aspect("equal")
    ax.set_axis_off()


def main(out):
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.6))
    panel(axes[0], THETA, f"(a) THETA +{THETA:g}: tilted away from the detector\n"
                          f"PSI = {THETA:g} + {PHI:g} = {THETA + PHI:g}°")
    panel(axes[1], -THETA, f"(b) THETA -{THETA:g}: tilted towards the detector\n"
                           f"PSI = |-{THETA:g} + {PHI:g}| = {abs(PHI - THETA):g}°")
    fig.text(
        0.5, 0.03,
        "GEOMETRY IBM, seen along the tilt axis: the axis is perpendicular to the scattering "
        "plane, so beam, detector and normal all stay in this plane\nand every angle is drawn "
        "true. PSI = |THETA + PHI|; the sign of THETA says which way the sample turns: "
        "away from the detector (+) or towards it (-).",
        ha="center", fontsize=10.5,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.1, wspace=0.05)
    fig.savefig(out, dpi=130)


if __name__ == "__main__":
    default = Path(__file__).resolve().parents[1] / "docs" / "assets" / "geometry-ibm.png"
    main(sys.argv[1] if len(sys.argv) > 1 else default)
