"""GEOMETRY GENERAL, in the Cornell figure's perspective: the red plane
stays horizontal, and the normal -- THETA from the beam -- rises RISE out of
it. Its dashed projection onto the red plane shows the swing, the dotted
line the rise. CORNELL is the case with no rise; in GENERAL PSI is typed in.

    PYTHONPATH=src python3.14 tools/geometry/general.py out.png
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
    BLUE, MAGENTA, ORANGE, PHI, RED, THETA,
    angle, arc2d, arc3d, arrow2d, arrow3d, rump_angles, unit,
)


RISE = 15.0  # the normal above the horizontal red plane

# x: the beam's travel; y: into the page; z: up -- as in the Cornell figure.
t, f, e = (math.radians(v) for v in (THETA, PHI, RISE))
SWING = math.degrees(math.acos(math.cos(t) / math.cos(e)))  # keeps THETA exact
w = math.radians(SWING)
S = np.array([-1.0, 0.0, 0.0])
D = np.array([-math.cos(f), 0.0, -math.sin(f)])
N = np.array([-math.cos(e) * math.cos(w), -math.cos(e) * math.sin(w), math.sin(e)])
FOOT = N * [1, 1, 0]                       # the normal's projection onto the red plane
PSI = angle(N, D)
CORNELL = math.degrees(math.acos(math.cos(t) * math.cos(f)))

def slab(ax, normal, half=(0.55, 0.85), thickness=0.12):
    u = unit(np.cross(normal, [0, 0, 1.0]))
    v = np.cross(u, normal)
    front = [a * half[0] * u + b * half[1] * v for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    back = [c - thickness * normal for c in front]
    faces = [front, back] + [[front[i], front[(i + 1) % 4], back[(i + 1) % 4], back[i]]
                             for i in range(4)]
    ax.add_collection3d(Poly3DCollection(
        faces, facecolors=["0.55", "0.35", "0.25", "0.3", "0.45", "0.3"],
        edgecolor="0.2", linewidths=0.6, alpha=0.5, zorder=1))


def overview(ax):
    ax.computed_zorder = False
    xs, zs = np.meshgrid([-2.4, 0.0], [-1.3, 0.6])
    ax.plot_surface(xs, np.zeros_like(xs), zs, color=ORANGE, alpha=0.16, lw=0, zorder=0)
    xs, ys = np.meshgrid([-2.4, 0.0], [-1.4, 0.7])
    ax.plot_surface(xs, ys, np.zeros_like(xs), color=RED, alpha=0.12, lw=0, zorder=0)
    slab(ax, N)
    arrow3d(ax, [-2.4, 0, 0], [2.33, 0, 0], "black")
    ax.text(-2.4, 0, 0.1, "beam", fontsize=10)
    arrow3d(ax, [0, 0, 0], 2.0 * D, MAGENTA)
    ax.text(*(2.05 * D + [0, 0, -0.15]), "detector", color=MAGENTA, fontsize=10,
            ha="center")
    length = 2.32
    arrow3d(ax, [0, 0, 0], length * N, BLUE)
    ax.text(*(0.8 * length * N + [0, 0, 0.12]), "normal", color=BLUE, fontsize=10,
            ha="center", va="bottom")
    # The normal's projection onto the red plane, and the rise above it.
    ax.plot(*np.array([[0, 0, 0], length * FOOT]).T, color=BLUE, ls="--", lw=1.6)
    ax.plot(*np.array([length * FOOT, length * N]).T, color=BLUE, ls=":", lw=1.4)
    arc3d(ax, N, S, 1.62, RED, "")
    arc3d(ax, S, D, 0.5 * 1.82, ORANGE, "")
    arc3d(ax, N, D, 1.15, MAGENTA, "")
    captions = [
        (RED, f"THETA {THETA:g}°: beam to normal"),
        (ORANGE, f"PHI {PHI:g}°: beam to detector, in the orange plane (beam and detector)"),
        (MAGENTA, f"PSI {PSI:.1f}°: normal to detector"),
        (BLUE, f"dashed: the normal's projection onto the red, horizontal plane --"),
        (BLUE, f"the normal rises {RISE:g}° out of it (dotted); in CORNELL it lies in it"),
    ]
    for i, (color, text) in enumerate(captions):
        ax.text2D(0.02, 0.17 - 0.045 * i, text, transform=ax.transAxes, color=color,
                  fontsize=9.5, fontweight="bold")
    ax.set_title(f"(a) Overview: PSI = {PSI:.1f}°", fontsize=11)
    ax.set_xlim(-2.8, 0.6)
    ax.set_ylim(-1.6, 1.0)
    ax.set_zlim(-1.3, 1.3)
    ax.set_box_aspect((3.0, 2.5, 2.6))
    ax.view_init(elev=25, azim=-105)
    ax.set_axis_off()


def flat(ax, title):
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-2.5, 1.0)
    ax.set_ylim(-1.7, 1.5)
    ax.set_aspect("equal")
    ax.set_axis_off()


def top_view(ax):
    """Looking down onto the red, horizontal plane: the swing is true size."""
    ax.add_patch(Polygon([(-2.45, -1.65), (0.95, -1.65), (0.95, 1.2), (-2.45, 1.2)],
                         facecolor=RED, alpha=0.07, edgecolor="none"))
    ax.text(-2.4, 1.1, "the red, horizontal plane", color=RED, fontsize=9, va="top")
    ax.text(-0.75, -1.6, "front (where you stand in (a) and (c))", color="0.45", fontsize=8,
            ha="center")

    def p(v):
        return np.array([v[0], v[1]])      # down the page = towards the front

    ax.plot([-2.4, 0.6], [0, 0], color=ORANGE, lw=5, alpha=0.35, solid_capstyle="butt")
    ax.text(-0.95, 0.06, "orange plane, edge-on", color=ORANGE, fontsize=8, va="bottom")
    arrow2d(ax, (-0.04, 0.0), "black", start=(-2.4, 0.0))
    ax.text(-2.4, 0.08, "beam", fontsize=10)
    arrow2d(ax, 1.6 * p(FOOT) / np.linalg.norm(p(FOOT)), BLUE, lw=1.6, ls="--")
    ax.text(*(1.65 * p(FOOT) / np.linalg.norm(p(FOOT)) + [0, -0.1]),
            f"normal, seen from above\n(it rises {RISE:g}° towards you)", color=BLUE,
            fontsize=9, ha="center", va="top")
    arrow2d(ax, 1.9 * p(D), "darkmagenta", lw=1.5, ls="--")
    ax.text(1.9 * p(D)[0], 0.3, "detector (below the beam,\nforeshortened)",
            color=MAGENTA, fontsize=8, va="bottom")
    arc2d(ax, p(FOOT), p(S), 0.9, "0.3", f"swing {SWING:.1f}°", at=(-1.05, -0.2), ha="right")
    flat(ax, "(b) Top view, onto the red plane:\nthe swing is true size")


def front_view(ax):
    """Looking onto the orange plane: PHI is true size."""
    ax.add_patch(Polygon([(-2.45, -1.65), (0.95, -1.65), (0.95, 1.45), (-2.45, 1.45)],
                         facecolor=ORANGE, alpha=0.08, edgecolor="none"))
    ax.text(-2.4, 1.35, "the plane of beam and detector (orange)", color=ORANGE, fontsize=9,
            va="top")

    def p(v):
        return np.array([v[0], v[2]])

    ax.plot([-2.4, 0.6], [0, 0], color=RED, lw=5, alpha=0.25, solid_capstyle="butt")
    ax.text(0.6, 0.08, "red plane,\nedge-on", color=RED, fontsize=8, ha="right", va="bottom")
    arrow2d(ax, (-0.04, 0.0), "black", start=(-2.4, 0.0))
    ax.text(-2.4, 0.08, "beam", fontsize=10)
    arrow2d(ax, 1.9 * p(D), MAGENTA)
    ax.text(*(1.95 * p(D) + [0, -0.12]), "detector", color=MAGENTA, fontsize=10,
            ha="center", va="top")
    arrow2d(ax, 1.6 * p(N), BLUE, lw=1.5, ls="--")
    ax.text(*(1.6 * p(N) + [0, 0.1]), "normal (rises above the red plane,\n"
            "leans towards you: foreshortened)", color=BLUE, fontsize=8, va="bottom",
            ha="center")
    arc2d(ax, p(D), p(S), 1.15, ORANGE, f"PHI {PHI:g}°", at=(-1.45, -0.42), ha="right")
    flat(ax, "(c) Front view, onto the orange plane:\nPHI is true size")


def commands(ax):
    inbound, outbound = rump_angles("GENERAL", THETA, round(PSI, 1))
    lines = [
        "RUMP:",
        "  GEOMETRY GENERAL",
        f"  THETA {THETA:g}     beam to the normal",
        f"  PHI {PHI:g}        180° - scattering angle ({180 - PHI:g}°)",
        f"  PSI {PSI:.1f}      detector to the normal: typed in",
        "",
        f"  -> pyRUMP uses beam {inbound:.0f}°, detector {outbound:.1f}°",
        "",
        "  For comparison, from THETA and PHI alone:",
        f"     CORNELL  PSI = {CORNELL:.1f}°   (no rise)",
        f"     IBM      PSI = {abs(PHI - THETA):g}° or {PHI + THETA:g}°",
    ]
    ax.text(0.05, 0.5, "\n".join(lines), transform=ax.transAxes, family="monospace",
            fontsize=11.5, va="center",
            bbox={"boxstyle": "round,pad=0.8", "facecolor": "white", "edgecolor": "0.6"})
    ax.set_axis_off()


def figure(out):
    fig = plt.figure(figsize=(14, 12.5))
    overview(fig.add_subplot(2, 2, 1, projection="3d"))
    top_view(fig.add_subplot(2, 2, 2))
    front_view(fig.add_subplot(2, 2, 3))
    commands(fig.add_subplot(2, 2, 4))
    fig.text(
        0.5, 0.02,
        f"GEOMETRY GENERAL: the normal can lie anywhere THETA from the beam. Here it rises "
        f"{RISE:g}° out of the horizontal plane where CORNELL keeps it,\nso PSI = {PSI:.1f}° "
        f"-- between CORNELL's {CORNELL:.1f}° and IBM's {PHI + THETA:g}° -- and no rule "
        "computes it: PSI is typed in.",
        ha="center", fontsize=10.5,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.07, wspace=0.05, hspace=0.12)
    fig.savefig(out, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    figure(sys.argv[1])
