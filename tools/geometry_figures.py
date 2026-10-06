"""Draw the figures of the manual's Experimental geometry page.

    PYTHONPATH=src python3.14 tools/geometry_figures.py

writes, into docs/assets/:

- geometry-ibm.png: GEOMETRY IBM seen along the tilt axis -- beam, detector
  and normal all lie in the scattering plane, so a flat drawing shows every
  angle true.
- geometry-cornell.png: GEOMETRY CORNELL as in a chamber -- beam horizontal,
  scattering plane and tilt axis vertical, the sample turning like a door: a
  3D overview with the beam-normal plane (red) and the beam-detector plane
  (orange), the top view (THETA true size), the front view (PHI true size)
  and the RUMP commands.

The exit angles in the command boxes are computed by pyRUMP itself
(``pyrump.model.geometry.Geometry``), so the figures can't drift from the
code. Both are drawn at THETA 25, PHI 35 -- larger than real RBS angles, so
the arcs stay apart.
"""

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Arc, Polygon  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402


# -- IBM: one flat view --------------------------------------------------------

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


def ibm_figure(out):
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


# -- 3D helpers ---------------------------------------------------------------


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def angle(a, b):
    return np.degrees(np.arccos(np.clip(np.dot(unit(a), unit(b)), -1, 1)))


def arc3d(ax, a, b, radius, color, label, shift=(0, 0, 0)):
    a, b = unit(a), unit(b)
    t = np.linspace(0, 1, 40)
    omega = np.arccos(np.clip(np.dot(a, b), -1, 1))
    pts = np.array([
        (np.sin((1 - s) * omega) * a + np.sin(s * omega) * b) / np.sin(omega) for s in t
    ]) * radius
    ax.plot(*pts.T, color=color, lw=2)
    mid = pts[len(pts) // 2] + np.asarray(shift, float)
    ax.text(*mid, label, color=color, fontsize=10, fontweight="bold", ha="center")


def arrow3d(ax, start, vec, color, lw=2.5, ls="-"):
    ax.quiver(*start, *vec, color=color, lw=lw, arrow_length_ratio=0.12, linestyle=ls)


# -- CORNELL: a 3D overview, the top and front views, the commands ------------

def arrow2d(ax, end, color, lw=2.5, ls="-", start=(0, 0)):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": lw, "ls": ls,
                            "mutation_scale": 16})


def arc2d(ax, a, b, radius, color, label, at, ha="center"):
    """The short arc between ``a`` and ``b`` (never the way round)."""
    t1, t2 = (math.degrees(math.atan2(v[1], v[0])) for v in (a, b))
    if t2 < t1:
        t1, t2 = t2, t1
    if t2 - t1 > 180:
        t1, t2 = t2, t1 + 360
    ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=t1, theta2=t2, color=color, lw=2))
    ax.text(*at, label, color=color, fontsize=11, fontweight="bold", ha=ha, va="center")


# x: the beam's direction of travel (left to right); y: into the page;
# z: up. The scattering plane is the vertical x-z plane (the page); the
# tilt axis is z, vertical, in that plane and across the beam.
t, f = math.radians(THETA), math.radians(PHI)
S = np.array([-1.0, 0.0, 0.0])                     # back towards the source
D = np.array([-math.cos(f), 0.0, -math.sin(f)])    # the RBS detector: back and down
N = np.array([-math.cos(t), -math.sin(t), 0.0])    # turned about z, towards the viewer
AXIS = np.array([0.0, 0.0, 1.0])
PSI = angle(N, D)

BLUE, MAGENTA, RED, ORANGE, PURPLE = "royalblue", "darkmagenta", "crimson", "darkorange", "purple"


def rump_psi():
    from pyrump.model.geometry import Geometry, GeometryKind

    g = Geometry(theta=THETA, phi=PHI, kind=GeometryKind.CORNELL)
    return math.degrees(math.acos(1.0 / g.sec_out))


# -- the overview --------------------------------------------------------------


def slab(ax, normal, half_width=0.55, half_height=0.85, thickness=0.12):
    """The sample: a plate whose front face, through the origin, faces
    ``normal``; the material lies behind it."""
    across = unit(np.cross(AXIS, normal))         # horizontal, in the face
    front = [a * half_width * across + b * half_height * AXIS
             for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    back = [c - thickness * normal for c in front]
    faces = [front, back] + [[front[i], front[(i + 1) % 4], back[(i + 1) % 4], back[i]]
                             for i in range(4)]
    shades = ["0.55", "0.35", "0.25", "0.3", "0.45", "0.3"]
    ax.add_collection3d(Poly3DCollection(faces, facecolors=shades, edgecolor="0.2",
                                         linewidths=0.6, alpha=0.5, zorder=1))


def rotation_mark(ax, height=1.15, radius=0.22):
    """An ellipse round the top of the tilt axis, with an arrowhead."""
    a = np.linspace(0.15 * np.pi, 1.85 * np.pi, 60)
    ring = np.array([radius * np.cos(a), radius * np.sin(a), np.full_like(a, height)])
    ax.plot(*ring, color=PURPLE, lw=1.5)
    end, before = ring[:, -1], ring[:, -4]
    ax.quiver(*before, *(end - before), color=PURPLE, lw=1.5, arrow_length_ratio=1.0)


def overview(ax):
    # Arrows and arcs are always drawn over the sample: matplotlib's own
    # depth sorting would hide them behind it.
    ax.computed_zorder = False
    # The beam-detector plane (the vertical scattering plane, where PHI
    # lies) and the beam-normal plane (horizontal here, where THETA lies).
    xs, zs = np.meshgrid([-2.4, 0.0], [-1.3, 0.6])
    ax.plot_surface(xs, np.zeros_like(xs), zs, color=ORANGE, alpha=0.16, lw=0, zorder=0)
    xs, ys = np.meshgrid([-2.4, 0.0], [-1.4, 0.7])
    ax.plot_surface(xs, ys, np.zeros_like(xs), color=RED, alpha=0.12, lw=0, zorder=0)

    slab(ax, N)
    ax.plot(*np.array([[0, 0, -1.3], [0, 0, 1.25]]).T, color=PURPLE, ls=":", lw=1.6)
    rotation_mark(ax)
    ax.text(0.08, 0, 1.3, "tilt axis", color=PURPLE, fontsize=9)

    arrow3d(ax, [-2.4, 0, 0], [2.33, 0, 0], "black")
    ax.text(-2.4, 0, 0.1, "beam", fontsize=10)
    arrow3d(ax, [0, 0, 0], 2.0 * D, MAGENTA)
    ax.text(*(2.05 * D + [0, 0, -0.15]), "RBS detector", color=MAGENTA, fontsize=10,
            ha="center")
    length = 0.8 * 2.9
    arrow3d(ax, [0, 0, 0], length * N, BLUE)
    ax.text(*(3.1 * N + [-0.15, -0.15, 0.15]), "normal", color=BLUE, fontsize=10, ha="right")

    arc3d(ax, N, S, 0.85 * 2 * 0.95, RED, "")
    arc3d(ax, S, D, 1.3 * 1.4, ORANGE, "")
    arc3d(ax, N, D, 1.15, MAGENTA, "")
    captions = [
        (RED, f"THETA {THETA:g}°: beam to normal, in the red plane (beam and normal)"),
        (ORANGE, f"PHI {PHI:g}°: beam to detector, in the orange plane (beam and detector)"),
        (MAGENTA, f"PSI {PSI:.1f}°: normal to detector, across the two planes"),
    ]
    for i, (color, text) in enumerate(captions):
        ax.text2D(0.02, 0.13 - 0.05 * i, text, transform=ax.transAxes, color=color,
                  fontsize=9.5, fontweight="bold")

    ax.set_title(f"(a) Overview: PSI = {PSI:.1f}°, true size in no flat view", fontsize=11)
    ax.set_xlim(-2.8, 0.6)
    ax.set_ylim(-1.6, 1.0)
    ax.set_zlim(-1.3, 1.3)
    ax.set_box_aspect((3.0, 2.5, 2.6))
    ax.view_init(elev=20, azim=-125)
    ax.set_axis_off()


# -- the flat views -------------------------------------------------------------


def flat(ax, title):
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-2.5, 1.0)
    ax.set_ylim(-1.7, 1.5)
    ax.set_aspect("equal")
    ax.set_axis_off()


def beam2d(ax):
    arrow2d(ax, (-0.04, 0.0), "black", start=(-2.4, 0.0))
    ax.text(-2.4, 0.08, "beam", fontsize=10)


def top_view(ax):
    """Looking down the vertical tilt axis: THETA is true size."""
    ax.add_patch(Polygon([(-2.45, -1.65), (0.95, -1.65), (0.95, 1.2), (-2.45, 1.2)],
                         facecolor=RED, alpha=0.07, edgecolor="none"))
    ax.text(-2.4, 1.1, "the plane of beam and normal (red)", color=RED, fontsize=9, va="top")
    ax.text(-0.75, -1.6, "front (where you stand in (a) and (c))", color="0.45", fontsize=8,
            ha="center")

    def p(v):
        return np.array([v[0], v[1]])      # down the page = towards the front

    across = unit(np.cross(AXIS, N))
    edge = [p(a * 0.55 * across + b * 0.12 * -N) for a, b in ((-1, 0), (1, 0), (1, 1), (-1, 1))]
    ax.add_patch(Polygon(edge, closed=True, facecolor="0.45", edgecolor="0.2"))
    ax.text(*(p(0.6 * across) + [0.0, -0.25]), "sample, seen from above:\nturned by THETA",
            color="0.3", fontsize=9, va="top", ha="center")
    ax.plot([-2.4, 0.6], [0, 0], color=ORANGE, lw=5, alpha=0.35, solid_capstyle="butt")
    ax.text(-0.95, 0.06, "orange plane, edge-on", color=ORANGE, fontsize=8, va="bottom")
    beam2d(ax)
    arrow2d(ax, 1.5 * p(N), BLUE)
    ax.text(*(1.58 * p(N) + [0, -0.08]), "normal", color=BLUE, fontsize=10, ha="right",
            va="top")
    arrow2d(ax, 1.9 * p(D), MAGENTA, lw=1.5, ls="--")
    ax.text(1.9 * p(D)[0], 0.3, "detector (below the beam,\nforeshortened)",
            color=MAGENTA, fontsize=8, va="bottom")
    arc2d(ax, p(N), p(S), 0.8, RED, f"THETA {THETA:g}°", at=(-1.0, -0.22), ha="right")
    ax.plot(0, 0, marker="o", ms=10, mfc="white", mec=PURPLE, mew=2)
    ax.plot(0, 0, marker=".", ms=5, color=PURPLE)
    ax.text(0.2, 0.12, "tilt axis\n(vertical, towards you)", color=PURPLE, fontsize=9,
            va="bottom")
    flat(ax, "(b) Top view, down the tilt axis:\nTHETA is true size")


def front_view(ax):
    """Looking across the vertical scattering plane: PHI is true size."""
    ax.add_patch(Polygon([(-2.45, -1.65), (0.95, -1.65), (0.95, 1.45), (-2.45, 1.45)],
                         facecolor=ORANGE, alpha=0.08, edgecolor="none"))
    ax.text(-2.4, 1.35, "the plane of beam and detector (orange)", color=ORANGE, fontsize=9,
            va="top")

    def p(v):
        return np.array([v[0], v[2]])

    across = unit(np.cross(AXIS, N))
    face = [p(a * 0.55 * across + b * 0.85 * AXIS) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ax.add_patch(Polygon(face, closed=True, facecolor="0.8", edgecolor="0.4", alpha=0.85))
    right = max(c[0] for c in face)
    ax.text(right + 0.25, 0.8, "sample,\nturned towards\nyou", color="0.35", fontsize=8,
            va="top")
    beam2d(ax)
    arrow2d(ax, 1.9 * p(D), MAGENTA)
    ax.text(*(1.95 * p(D) + [0, -0.12]), "RBS detector", color=MAGENTA, fontsize=10,
            ha="center", va="top")
    arrow2d(ax, 1.5 * p(N), BLUE, lw=1.5, ls="--")
    ax.text(1.5 * p(N)[0], 0.12, "normal (turned towards you:\nforeshortened)",
            color=BLUE, fontsize=8, va="bottom")
    arc2d(ax, p(D), p(S), 1.15, ORANGE, f"PHI {PHI:g}°", at=(-1.45, -0.42), ha="right")
    ax.plot([0, 0], [-1.3, 1.25], color=PURPLE, ls=":", lw=1.6)
    ax.text(0.06, 1.2, "tilt axis", color=PURPLE, fontsize=9)
    flat(ax, "(c) Front view, across the scattering plane:\nPHI is true size")


def commands(ax):
    lines = [
        "RUMP:",
        "  GEOMETRY CORNELL",
        f"  THETA {THETA:g}     beam to the normal",
        f"  PHI {PHI:g}        180° - scattering angle ({180 - PHI:g}°)",
        "  (PSI is not used: RUMP computes it)",
        "",
        "  -> detector to the normal:",
        f"     PSI = acos(cos THETA x cos PHI) = {rump_psi():.1f}°",
    ]
    ax.text(0.08, 0.5, "\n".join(lines), transform=ax.transAxes, family="monospace",
            fontsize=12, va="center",
            bbox={"boxstyle": "round,pad=0.8", "facecolor": "white", "edgecolor": "0.6"})
    ax.set_axis_off()


def cornell_figure(out):
    fig = plt.figure(figsize=(14, 12.5))
    overview(fig.add_subplot(2, 2, 1, projection="3d"))
    top_view(fig.add_subplot(2, 2, 2))
    front_view(fig.add_subplot(2, 2, 3))
    commands(fig.add_subplot(2, 2, 4))
    fig.text(
        0.5, 0.02,
        "GEOMETRY CORNELL: the tilt axis lies in the scattering plane, across the beam -- here "
        "both vertical -- so the sample turns like a door\nand its normal swings out of the "
        "scattering plane. THETA -25 turns it the other way: the mirror image, the same PSI. "
        "A typed PSI is ignored.",
        ha="center", fontsize=10.5,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.07, wspace=0.05, hspace=0.12)
    fig.savefig(out, dpi=130)


def main(folder: Path) -> None:
    ibm_figure(folder / "geometry-ibm.png")
    cornell_figure(folder / "geometry-cornell.png")


if __name__ == "__main__":
    default = Path(__file__).resolve().parents[1] / "docs" / "assets"
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else default)
