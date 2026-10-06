"""What the geometry figures share: the angles they are drawn at, the colour
code, and the arrows and arcs, flat and 3D.

The figures say "detector", not "RBS detector": the same three angles --
beam to normal (THETA), beam to detector (PHI), normal to detector (PSI) --
describe the RBS detector and the PIXE detector alike.
"""

import math

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
from matplotlib.patches import Arc  # noqa: E402

#: Drawn larger than real RBS angles, so the arcs stay apart.
THETA, PHI = 25.0, 35.0

BLUE, MAGENTA, RED, ORANGE, PURPLE = "royalblue", "darkmagenta", "crimson", "darkorange", "purple"


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def angle(a, b) -> float:
    """The angle between two directions, in degrees."""
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(a), unit(b)), -1, 1))))


def rump_angles(kind: str, theta: float, psi: float = 0.0) -> tuple[float, float]:
    """The beam and detector angles to the normal that pyRUMP itself takes
    from ``GEOMETRY kind``, ``THETA theta``, ``PHI`` and ``PSI psi`` -- so a
    figure's command box can't drift from the code."""
    from pyrump.model.geometry import Geometry, GeometryKind

    geometry = Geometry(theta=theta, phi=PHI, psi=psi, kind=GeometryKind[kind])
    return (math.degrees(math.acos(1.0 / geometry.sec_in)),
            math.degrees(math.acos(1.0 / geometry.sec_out)))


# -- flat drawings -------------------------------------------------------------


def arrow2d(ax, end, color, lw=2.5, ls="-", start=(0, 0)):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": lw, "ls": ls,
                            "mutation_scale": 16})


def arc2d(ax, a, b, radius, color, label, at, ha="center"):
    """The short arc between directions ``a`` and ``b`` (never the way
    round), labelled at the point ``at``."""
    t1, t2 = (math.degrees(math.atan2(v[1], v[0])) for v in (a, b))
    if t2 < t1:
        t1, t2 = t2, t1
    if t2 - t1 > 180:
        t1, t2 = t2, t1 + 360
    ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=t1, theta2=t2, color=color, lw=2))
    ax.text(*at, label, color=color, fontsize=11, fontweight="bold", ha=ha, va="center")


# -- 3D drawings ---------------------------------------------------------------


def arc3d(ax, a, b, radius, color, label="", shift=(0, 0, 0)):
    """The arc between directions ``a`` and ``b`` on a sphere of ``radius``."""
    a, b = unit(a), unit(b)
    omega = np.arccos(np.clip(np.dot(a, b), -1, 1))
    pts = np.array([
        (np.sin((1 - s) * omega) * a + np.sin(s * omega) * b) / np.sin(omega)
        for s in np.linspace(0, 1, 40)
    ]) * radius
    ax.plot(*pts.T, color=color, lw=2)
    if label:
        mid = pts[len(pts) // 2] + np.asarray(shift, float)
        ax.text(*mid, label, color=color, fontsize=10, fontweight="bold", ha="center")


def arrow3d(ax, start, vec, color, lw=2.5, ls="-"):
    ax.quiver(*start, *vec, color=color, lw=lw, arrow_length_ratio=0.12, linestyle=ls)
