## Experimental geometry

An RBS spectrum depends on three directions: the **beam**, the **detector**
and the **sample normal**. Three angles between them describe the
experiment, each set on the data buffer:

| Command | Angle between | What it sets |
|---|---|---|
| [`THETA`](buffers.md#theta) | the beam (looking back at the source) and the sample normal | the beam's path into the sample |
| [`PHI`](buffers.md#phi) | the beam (looking back at the source) and the detector: **180° minus the scattering angle** | the kinematics and the cross section |
| [`PSI`](buffers.md#psi) | the detector and the sample normal | the particles' path out of the sample |

A detector at a scattering angle of 170° is `PHI 10`, not `PHI 170` — RUMP's
convention, kept by pyRUMP (see [Things that will catch you
out](gotchas.md)).

The three angles are not independent: they are the angles between three
directions, so `PSI` can only lie between |THETA − PHI| and THETA + PHI.
[`GEOMETRY`](buffers.md#geometry) says how pyRUMP gets `PSI`, and there are
three cases:

| `GEOMETRY` | Where the detector is | `PSI` |
|---|---|---|
| [`IBM`](#ibm-geometry) | beside the beam, in the plane the sample normal turns in | computed: PSI = \|THETA + PHI\| |
| [`CORNELL`](#cornell-geometry) | above or below the beam, in the plane of the tilt axis | computed: cos PSI = cos THETA × cos PHI |
| [`GENERAL`](#general-geometry) | anywhere | typed in |

With `IBM` or `CORNELL`, a typed `PSI` is ignored.

The figures below show all three in the same chamber, at the same angles
(`THETA 25`, `PHI 35`), and with the same colours: the beam in black, the
detector in purple, the sample normal in blue; `THETA` in red, `PHI` in
orange, `PSI` in purple. In each, the beam comes in horizontally and the
sample turns about a vertical tilt axis, like a door: only the detector
moves from one case to the next.

### IBM geometry

The detector sits beside the beam, in the horizontal plane (red) in which
the sample normal turns. The tilt axis is then perpendicular to the
scattering plane — the plane of the beam and the detector — and the beam,
the detector and the normal never leave the red plane. Looking down the tilt
axis, every angle is drawn true size. The figure shows the two ways the
sample can turn, each with its 3D view (a, b), its top view (c, d) and its
RUMP commands (e, f):

![GEOMETRY IBM at PHI 35, the sample turning about a vertical tilt axis with the detector beside the beam: (a) THETA +25, turned away from the detector, PSI = 60 degrees; (b) THETA -25, turned towards it, PSI = 10 degrees; (c, d) their top views, every angle true size; (e, f) their RUMP commands](../assets/geometry-ibm.png)

PSI = |THETA + PHI|. The sign of `THETA` says which way the sample turns:
positive away from the detector, negative towards it. RUMP's own manual:
"The sign ambiguity is resolved so that PSI = THETA+PHI. If the sample normal
points to the detector, then THETA should equal −PHI."

`examples/MnPt.RBS` has `GEOMETRY GENERAL`, `THETA -9`, `PHI 11`, `PSI 20`.
Its PSI 20 = 9 + 11 is the IBM value for a sample turned 9° away from the
detector, which RUMP's IBM convention would write as `THETA 9`. The
acquisition software writes `THETA -9` and spells `PSI` out, so the file
needs `GENERAL`.

### Cornell geometry

The detector sits below the beam, so the tilt axis lies in the scattering
plane, across the beam. In the figure both are vertical: the beam comes in
horizontally, the detector looks back and down, and the sample turns about
the vertical axis like a door. Its normal therefore swings *out of* the
scattering plane, and the three angles lie in different planes:

* `THETA` in the plane of the beam and the normal (red, horizontal here);
* `PHI` in the plane of the beam and the detector, the scattering plane
  (orange, vertical);
* `PSI` across the two, between the normal and the detector.

No single flat drawing shows all three at their true size. The overview (a)
shows them together; the top view (b), looking down the tilt axis onto the
red plane, shows `THETA` true; the front view (c), looking onto the orange
plane, shows `PHI` true:

![GEOMETRY CORNELL at THETA 25, PHI 35, the sample turning about a vertical tilt axis: (a) a 3D overview with the beam-normal plane in red, the beam-detector plane in orange and PSI = 42.1 degrees; (b) the top view, THETA 25 degrees true size; (c) the front view, PHI 35 degrees true size; and the RUMP commands](../assets/geometry-cornell.png)

cos PSI = cos THETA × cos PHI. The sign of `THETA` only says which way the
sample turns, so `THETA -25` gives the same `PSI` as `THETA 25`. `PSI`
always lies between the two IBM values for the same `THETA` and `PHI` —
here 42.1°, between 10° and 60°.

### General geometry

With `GENERAL`, pyRUMP uses `PSI` as typed: it is for a normal that neither
rule describes — a sample tilted about two axes, or a detector out of both
the IBM and the CORNELL positions. In the figure the normal is still
`THETA` 25° from the beam, but it rises 15° out of the red plane where
CORNELL keeps it. Its dashed projection onto the red plane shows how far the
sample has swung (20.2°, true size in the top view b); the dotted line, how
far it rises. `PSI` comes out at 53.6° — between CORNELL's 42.1° and IBM's
60° — and only a typed `PSI` can say so:

![GEOMETRY GENERAL at THETA 25, PHI 35, with the normal rising 15 degrees out of the horizontal red plane: (a) a 3D overview, PSI = 53.6 degrees; (b) the top view, the 20.2-degree swing true size; (c) the front view, PHI 35 degrees true size; and the RUMP commands with PSI typed in](../assets/geometry-general.png)

Acquisition files often use `GENERAL` and spell all three angles out:
`examples/MnPt.RBS` has `GEOMETRY GENERAL`, `THETA -9`, `PHI 11`, `PSI 20`
(see [IBM geometry](#ibm-geometry) for what its `PSI` means).

### Fitting the tilt

When fitting, [`PERT THETA`](pert.md#theta) varies `THETA`; with `IBM` or
`CORNELL` the exit angle follows it, with `GENERAL` it stays at `PSI`.
