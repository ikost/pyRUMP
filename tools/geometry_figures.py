"""Draw the figures of the manual's Experimental geometry page.

    PYTHONPATH=src python3.14 tools/geometry_figures.py

writes, into docs/assets/:

- geometry-ibm.png: GEOMETRY IBM seen along the tilt axis -- beam, detector
  and normal all lie in the scattering plane, so a flat drawing shows every
  angle true.
- geometry-cornell.png: GEOMETRY CORNELL, whose normal leaves the scattering
  plane: a 3D overview, the views in which THETA and PHI are true size, and
  the RUMP commands.

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


def plate(ax, normal, size=0.75):
    n = unit(normal)
    u = unit(np.cross(n, [0, 1, 0]) if abs(n[1]) < 0.9 else np.cross(n, [1, 0, 0]))
    v = np.cross(n, u)
    corners = [s * u * size + t * v * size for s, t in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ax.add_collection3d(Poly3DCollection([corners], facecolor="0.55", alpha=0.45, edgecolor="0.3"))


# -- CORNELL: a 3D overview and two flat views --------------------------------

# x: across the beam, in the scattering plane (the tilt axis); y: out of the
# scattering plane; z: back towards the beam's source.
S = np.array([0, 0, 1.0])
D = unit([-math.sin(math.radians(PHI)), 0, math.cos(math.radians(PHI))])
N = np.array([0, -math.sin(math.radians(THETA)), math.cos(math.radians(THETA))])
PSI = angle(N, D)

BLUE, MAGENTA, RED, ORANGE, PURPLE = "royalblue", "darkmagenta", "crimson", "darkorange", "purple"


def rump_psi():
    from pyrump.model.geometry import Geometry, GeometryKind

    g = Geometry(theta=THETA, phi=PHI, kind=GeometryKind.CORNELL)
    return math.degrees(math.acos(1.0 / g.sec_out))


# -- flat views ---------------------------------------------------------------


def arrow2d(ax, end, color, lw=2.5, ls="-", start=(0, 0)):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": lw, "ls": ls,
                            "mutation_scale": 16})


def arc2d(ax, a, b, radius, color, label, at, ha="center"):
    p1, p2 = (math.degrees(math.atan2(v[1], v[0])) for v in (a, b))
    t1, t2 = sorted((p1, p2))
    ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=t1, theta2=t2, color=color, lw=2))
    ax.text(*at, label, color=color, fontsize=11, fontweight="bold", ha=ha, va="center")


def flat(ax, title):
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.3, 2.3)
    ax.set_aspect("equal")
    ax.set_axis_off()


def sample_outline(ax, project):
    """The sample plate, projected."""
    u = np.array([1.0, 0, 0])               # the tilt axis stays in the plate
    v = np.cross(N, u)
    corners = [project(a * 0.95 * u + b * 0.95 * v) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ax.add_patch(Polygon(corners, closed=True, facecolor="0.8", edgecolor="0.4", alpha=0.8))


def view_along_tilt_axis(ax):
    """Looking along x: the plane of beam and normal. THETA is true here."""
    def p(v):
        return np.array([-v[1], v[2]])     # right = out of the scattering plane

    sample_outline(ax, p)
    arrow2d(ax, (0, 0.04), "black", start=(0, 2.1))
    ax.text(0.06, 2.0, "beam", fontsize=10)
    arrow2d(ax, 1.4 * p(N), BLUE)
    ax.text(*(1.48 * p(N)), "normal", color=BLUE, fontsize=10)
    arrow2d(ax, 1.75 * p(D), MAGENTA, lw=1.5, ls="--")
    ax.text(-0.08, 1.75 * p(D)[1], "detector\n(behind the beam,\nforeshortened)",
            color=MAGENTA, fontsize=8, ha="right", va="center")
    arc2d(ax, p(N), p(S), 0.6, RED, f"THETA {THETA:g}°", at=(0.3, 0.75), ha="left")
    ax.plot(0, 0, marker="o", ms=10, mfc="white", mec=PURPLE, mew=2)
    ax.plot(0, 0, marker=".", ms=5, color=PURPLE)
    ax.text(0.1, -0.3, "tilt axis (towards you)", color=PURPLE, fontsize=9)
    edge = 0.95 * np.array([N[2], N[1]])
    ax.text(*(edge + [0.05, -0.12]), "sample, edge-on", color="0.35", fontsize=9, ha="left")
    flat(ax, "(b) Seen along the tilt axis:\nTHETA is true size")


def view_across_plane(ax):
    """Looking along y: the scattering plane. PHI is true here."""
    def p(v):
        return np.array([v[0], v[2]])

    sample_outline(ax, p)
    arrow2d(ax, (0, 0.04), "black", start=(0, 2.1))
    ax.text(0.06, 2.0, "beam", fontsize=10)
    arrow2d(ax, 1.75 * p(D), MAGENTA)
    ax.text(*(1.82 * p(D)), "RBS detector", color=MAGENTA, fontsize=10, ha="right")
    arrow2d(ax, 1.4 * p(N), BLUE, lw=1.5, ls="--")
    ax.text(0.06, 1.4 * p(N)[1] - 0.12, "normal (leans towards\nyou: foreshortened)",
            color=BLUE, fontsize=8, ha="left", va="top")
    arc2d(ax, p(S), p(D), 0.9, ORANGE, f"PHI {PHI:g}°", at=(-0.3, 1.08), ha="center")
    ax.plot([-1.25, 1.25], [0, 0], color=PURPLE, ls="-.", lw=1.5)
    ax.text(1.28, 0.0, "tilt axis", color=PURPLE, fontsize=9, va="center")
    flat(ax, "(c) Seen across the scattering plane:\nPHI is true size")


# -- the 3D overview -----------------------------------------------------------


def overview(ax):
    xs, zs = np.meshgrid([-1.6, 1.4], [-0.4, 2.2])
    ax.plot_surface(xs, np.zeros_like(xs), zs, color="orange", alpha=0.10, lw=0)
    ax.text(1.0, 0, 2.05, "scattering plane", color=ORANGE, fontsize=8)
    plate(ax, N)
    arrow3d(ax, [0, 0, 2.3], [0, 0, -2.1], "black")
    ax.text(0.08, 0, 2.3, "beam", fontsize=10)
    arrow3d(ax, [0, 0, 0], 1.95 * D, MAGENTA)
    ax.text(*(2.0 * D + [-0.1, 0, 0.08]), "RBS detector", color=MAGENTA, fontsize=10, ha="right")
    arrow3d(ax, [0, 0, 0], 1.4 * N, BLUE)
    ax.text(*(1.5 * N + [0, -0.2, 0.15]), "normal", color=BLUE, fontsize=10)
    tip = 1.4 * N
    ax.plot(*np.array([tip, tip * [1, 0, 1]]).T, color=BLUE, ls=":", lw=1.5)
    ax.plot(*np.array([[-1.15, 0, 0], [1.15, 0, 0]]).T, color=PURPLE, ls="-.", lw=1.8)
    ax.text(1.2, 0, 0, "tilt axis", color=PURPLE, fontsize=9)
    arc3d(ax, N, S, 0.6, RED, f"THETA {THETA:g}°", shift=(0.1, -0.75, 0.05))
    arc3d(ax, S, D, 1.05, ORANGE, f"PHI {PHI:g}°", shift=(-0.15, 0, 0.3))
    arc3d(ax, N, D, 0.85, MAGENTA, f"PSI {PSI:.1f}°", shift=(-0.55, 0.3, -0.35))
    ax.set_title(f"(a) Overview: PSI = {PSI:.1f}°, true size in no flat view", fontsize=11)
    ax.set_xlim(-1.6, 1.4)
    ax.set_ylim(-1.3, 1.3)
    ax.set_zlim(-0.6, 2.3)
    ax.set_box_aspect((3.0, 2.6, 2.9))
    ax.view_init(elev=22, azim=-110)
    ax.set_axis_off()


def commands(ax):
    """The RUMP commands for this geometry, where a fourth view would go."""
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
    view_along_tilt_axis(fig.add_subplot(2, 2, 2))
    view_across_plane(fig.add_subplot(2, 2, 3))
    commands(fig.add_subplot(2, 2, 4))
    lines = [
        "GEOMETRY CORNELL: the tilt axis lies in the scattering plane, across the beam, so the normal tilts out "
        "of the plane.\nTHETA -25 tilts it to the other side: "
        "the mirror image, and the same PSI. A typed PSI is ignored.",
    ]
    fig.text(0.5, 0.02, "\n".join(lines), ha="center", fontsize=10.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.07, wspace=0.05, hspace=0.12)
    fig.savefig(out, dpi=130)


def main(folder: Path) -> None:
    ibm_figure(folder / "geometry-ibm.png")
    cornell_figure(folder / "geometry-cornell.png")


if __name__ == "__main__":
    default = Path(__file__).resolve().parents[1] / "docs" / "assets"
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else default)
