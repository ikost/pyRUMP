"""GEOMETRY IBM in the Cornell figure's chamber: the beam horizontal, the
sample turning about a vertical tilt axis like a door -- as for CORNELL --
but the detector beside the beam, in the horizontal plane (red). The tilt
axis is then perpendicular to the scattering plane, and beam, detector and
normal all stay in the red plane. One column for each sign of THETA: the
3D overview, the top view (every angle true size) and the RUMP commands.

    PYTHONPATH=src python3.14 tools/geometry/ibm.py out.png
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from common import (  # noqa: E402
    BLUE, MAGENTA, ORANGE, PHI, PURPLE, RED, THETA,
    angle, arc2d, arc3d, arrow2d, arrow3d, rump_angles, unit,
)

# x: the beam's travel; y: into the page; z: up -- as in the Cornell figure.
# The tilt axis is z, vertical; the detector lies in the horizontal x-y
# plane, on the far side of the beam.
f = math.radians(PHI)
S = np.array([-1.0, 0.0, 0.0])                    # back towards the source
D = np.array([-math.cos(f), math.sin(f), 0.0])    # the detector: back and to the far side
AXIS = np.array([0.0, 0.0, 1.0])


def normal(theta: float) -> np.ndarray:
    """IBM's normal: turned by THETA about the vertical axis, in the red
    plane -- towards the viewer, away from the detector, for THETA > 0
    (PSI = THETA + PHI)."""
    t = math.radians(theta)
    return np.array([-math.cos(t), -math.sin(t), 0.0])


def slab(ax, n, half_width=0.55, half_height=0.85, thickness=0.12):
    """The sample, as in the Cornell figure: an upright plate facing ``n``."""
    across = unit(np.cross(AXIS, n))
    front = [a * half_width * across + b * half_height * AXIS
             for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    back = [c - thickness * n for c in front]
    faces = [front, back] + [[front[i], front[(i + 1) % 4], back[(i + 1) % 4], back[i]]
                             for i in range(4)]
    ax.add_collection3d(Poly3DCollection(
        faces, facecolors=["0.55", "0.35", "0.25", "0.3", "0.45", "0.3"],
        edgecolor="0.2", linewidths=0.6, alpha=0.5, zorder=1))


def rotation_mark(ax, height=1.15, radius=0.22):
    """An ellipse round the top of the tilt axis, with an arrowhead."""
    a = np.linspace(0.15 * np.pi, 1.85 * np.pi, 60)
    ring = np.array([radius * np.cos(a), radius * np.sin(a), np.full_like(a, height)])
    ax.plot(*ring, color=PURPLE, lw=1.5)
    end, before = ring[:, -1], ring[:, -4]
    ax.quiver(*before, *(end - before), color=PURPLE, lw=1.5, arrow_length_ratio=1.0)


def overview(ax, theta, title):
    ax.computed_zorder = False
    n = normal(theta)
    psi = angle(n, D)
    xs, ys = np.meshgrid([-2.4, 0.3], [-1.4, 1.4])
    ax.plot_surface(xs, ys, np.zeros_like(xs), color=RED, alpha=0.12, lw=0, zorder=0)
    slab(ax, n)
    ax.plot(*np.array([[0, 0, -1.3], [0, 0, 1.25]]).T, color=PURPLE, ls=":", lw=1.6)
    rotation_mark(ax)
    ax.text(0.08, 0, 1.3, "tilt axis", color=PURPLE, fontsize=9)
    arrow3d(ax, [-2.4, 0, 0], [2.33, 0, 0], "black")
    ax.text(-2.4, 0, 0.1, "beam", fontsize=10)
    arrow3d(ax, [0, 0, 0], 2.0 * D, MAGENTA)
    ax.text(*(2.1 * D + [0, 0, 0.08]), "detector", color=MAGENTA, fontsize=10, ha="center")
    length = 2.32
    arrow3d(ax, [0, 0, 0], length * n, BLUE)
    ax.text(*(1.12 * length * n), "normal", color=BLUE, fontsize=10, ha="center",
            va="top")
    arc3d(ax, n, S, 1.62, RED)
    arc3d(ax, S, D, 0.91, ORANGE)
    arc3d(ax, n, D, 1.15, MAGENTA)
    captions = [
        (RED, f"THETA {theta:+g}°: beam to normal"),
        (ORANGE, f"PHI {PHI:g}°: beam to detector"),
        (MAGENTA, f"PSI {psi:.0f}°: normal to detector"),
    ]
    for i, (color, text) in enumerate(captions):
        ax.text2D(0.02, 0.12 - 0.045 * i, text, transform=ax.transAxes, color=color,
                  fontsize=9.5, fontweight="bold")
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-2.8, 0.6)
    ax.set_ylim(-1.6, 1.4)
    ax.set_zlim(-1.3, 1.3)
    ax.set_box_aspect((3.0, 2.8, 2.6))
    ax.view_init(elev=40, azim=-100)
    ax.set_axis_off()


def top_view(ax, theta, title):
    """Looking down the vertical tilt axis, onto the red plane: every angle
    true size. The front -- where the viewer of (a) stands -- at the bottom."""
    ax.add_patch(Polygon([(-2.45, -1.65), (0.95, -1.65), (0.95, 1.55), (-2.45, 1.55)],
                         facecolor=RED, alpha=0.07, edgecolor="none"))
    ax.text(-2.4, 1.45, "the red, horizontal plane: beam, detector and normal", color=RED,
            fontsize=9, va="top")
    ax.text(-0.75, -1.6, "front (where you stand in (a) and (b))", color="0.45", fontsize=8,
            ha="center")

    def p(v):
        return np.array([v[0], v[1]])

    n = normal(theta)
    psi = angle(n, D)
    across = unit(np.cross(AXIS, n))
    edge = [p(a * 0.85 * across + b * 0.12 * -n) for a, b in ((-1, 0), (1, 0), (1, 1), (-1, 1))]
    ax.add_patch(Polygon(edge, closed=True, facecolor="0.55", edgecolor="0.25"))
    arrow2d(ax, (-0.04, 0.0), "black", start=(-2.4, 0.0))
    ax.text(-2.4, 0.08, "beam", fontsize=10)
    arrow2d(ax, 1.9 * p(D), MAGENTA)
    ax.text(*(1.95 * p(D) + [0, 0.1]), "detector", color=MAGENTA, fontsize=10, ha="center",
            va="bottom")
    arrow2d(ax, 1.6 * p(n), BLUE)
    ax.text(*(1.65 * p(n) + [0, -0.08 if theta > 0 else 0.08]), "normal", color=BLUE,
            fontsize=10, ha="center", va="top" if theta > 0 else "bottom")
    ax.plot(0, 0, marker="o", ms=10, mfc="white", mec=PURPLE, mew=2)
    ax.plot(0, 0, marker=".", ms=5, color=PURPLE)
    ax.text(0.15, 0.12 if theta > 0 else -0.12, "tilt axis\n(vertical, towards you)",
            color=PURPLE, fontsize=8, va="bottom" if theta > 0 else "top")
    # Three arcs meet at one vertex: distinct radii, labels in free space.
    if theta > 0:
        arc2d(ax, p(n), p(S), 0.55, RED, f"THETA {THETA:g}°", at=(-0.68, -0.06), ha="right")
        ax.texts[-1].set_va("top")
        arc2d(ax, p(S), p(D), 0.95, ORANGE, f"PHI {PHI:g}°", at=(-1.05, 0.33), ha="right")
        arc2d(ax, p(n), p(D), 1.35, MAGENTA, f"PSI {psi:.0f}°", at=(-1.49, -0.13), ha="right")
        ax.texts[-1].set_va("top")
    else:
        arc2d(ax, p(n), p(S), 0.75, RED, f"THETA {THETA:g}°", at=(-0.72, -0.06), ha="center")
        ax.texts[-1].set_va("top")
        arc2d(ax, p(S), p(D), 1.15, ORANGE, f"PHI {PHI:g}°", at=(-1.3, 0.38), ha="right")
        arc2d(ax, p(n), p(D), 0.5, MAGENTA, f"PSI {psi:.0f}°", at=(-0.38, 0.66), ha="center")
        ax.texts[-1].set_va("bottom")
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-2.5, 1.0)
    ax.set_ylim(-1.7, 1.6)
    ax.set_aspect("equal")
    ax.set_axis_off()


def commands(ax, theta, label):
    psi = rump_angles("IBM", theta)[1]
    way = "away from" if theta > 0 else "towards"
    lines = [
        "RUMP:",
        "  GEOMETRY IBM",
        f"  {f'PHI {PHI:g}':<12}180° - scattering angle ({180 - PHI:g}°)",
        f"  {f'THETA {theta:g}':<12}beam to the normal, turned",
        f"  {'':<12}{way} the detector",
        "  (PSI is not used: RUMP computes it)",
        "",
        f"  -> PSI = |THETA + PHI| = |{theta:g} + {PHI:g}| = {psi:.0f}°",
    ]
    ax.text(0.5, 0.5, "\n".join(lines), transform=ax.transAxes, family="monospace",
            fontsize=11.5, ha="center", va="center", multialignment="left",
            bbox={"boxstyle": "round,pad=0.8", "facecolor": "white", "edgecolor": "0.6"})
    ax.set_title(label, fontsize=11)
    ax.set_axis_off()


def figure(out):
    away = rump_angles("IBM", THETA)[1]
    towards = rump_angles("IBM", -THETA)[1]
    fig = plt.figure(figsize=(14, 16.5))
    grid = fig.add_gridspec(3, 2, height_ratios=(1, 1, 0.42))
    overview(fig.add_subplot(grid[0, 0], projection="3d"), THETA,
             f"(a) THETA +{THETA:g}: turned away from the detector\n"
             f"PSI = |{THETA:g} + {PHI:g}| = {away:.0f}°")
    overview(fig.add_subplot(grid[0, 1], projection="3d"), -THETA,
             f"(b) THETA -{THETA:g}: turned towards the detector\n"
             f"PSI = |-{THETA:g} + {PHI:g}| = {towards:.0f}°")
    top_view(fig.add_subplot(grid[1, 0]), THETA, f"(c) Top view of (a), THETA +{THETA:g}")
    top_view(fig.add_subplot(grid[1, 1]), -THETA, f"(d) Top view of (b), THETA -{THETA:g}")
    commands(fig.add_subplot(grid[2, 0]), THETA, f"(e) RUMP commands for (a)")
    commands(fig.add_subplot(grid[2, 1]), -THETA, f"(f) RUMP commands for (b)")
    fig.text(
        0.5, 0.012,
        "GEOMETRY IBM: the sample turns about a vertical axis as for CORNELL, but the detector sits beside "
        "the beam: the tilt axis is perpendicular to the scattering plane,\nso beam, detector "
        "and normal all stay in the red, horizontal plane, every angle is true size from above, "
        "and the sign of THETA says which way.",
        ha="center", fontsize=10.5,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.05, wspace=0.05, hspace=0.08)
    fig.savefig(out, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    figure(sys.argv[1])
