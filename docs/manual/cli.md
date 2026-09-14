## CLI reference

```
pyrump [--data DIR] {shell,simulate,fit,convert,plot} ...
```

### `pyrump simulate`

Simulate a spectrum from a sample description.

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
```

| Option | Default | Meaning |
|---|---|---|
| `--beam` | `4He` | Beam species: `4He`, `1H`, `He`, `28Si`… |
| `--energy` | `2.0` | Beam energy, MeV |
| `--theta` | `0.0` | Incidence angle from the sample normal, degrees |
| `--phi` | `10.0` | **180° minus the scattering angle** (10 → 170°) |
| `--psi` | `0.0` | Exit angle from the normal (GENERAL geometry only) |
| `--geometry` | `cornell` | `cornell`, `ibm` or `general` |
| `--kevch` | `5.0` | keV per channel |
| `--kev0` | `0.0` | keV at channel zero |
| `--channels` | `1024` | Number of channels |
| `--fwhm` | `15.0` | Detector resolution, keV |
| `--omega` | `1.0` | Detector solid angle, msr |
| `--charge` | `10.0` | Integrated charge, µC |
| `--current` | `0.0` | Beam current, nA — enables pile-up with `--tau` |
| `--tau` | `0.0` | Shaping time, µs |
| `-o, --output` | stdout | `.rbs` for binary, anything else for ASCII |
| `--two-column` | off | Write `channel value` instead of one column |

Output extension chooses the format. With no `-o`, energy/counts pairs go to
stdout and a summary to stderr.

### `pyrump fit`

Adjust sample parameters until the simulation matches a measurement.

```bash
pyrump fit start.lcm measured.rbs --vary thickness:0 --window 190 226
```

Takes every `simulate` option, plus:

| Option | Meaning |
|---|---|
| `--vary NAME` | Parameter to fit; repeat for several |
| `--window LOW HIGH` | Channel range to fit over; repeat for up to 10 |

`--vary` accepts `thickness:N` for layer *N*, or any of `mev`, `theta`, `phi`,
`psi`, `fwhm`, `tau`, `current`, `correction`, `kev/ch`, `kev(0)`, `straggle`,
`multiple_scatter`.

**Choose the window deliberately.** It should cover the part of the spectrum
that constrains what you are fitting, and no more — see
[Things that will catch you out](gotchas.md).

### `pyrump convert`

```bash
pyrump convert measured.rbs measured.dat        # binary → ASCII
pyrump convert spectrum.dat spectrum.txt --two-column
```

Writing `.rbs` requires beam and geometry metadata, so it only works from a
`.rbs` source — ASCII files do not carry it.

### `pyrump plot`

```bash
pyrump plot measured.rbs                            # interactive window
pyrump plot measured.rbs -o spectrum.png            # save
pyrump plot measured.rbs --compare simulated.rbs -o comparison.png
```

With `--compare` you get the data, the simulation over it, and a residual panel
showing the **Poisson residuals the fit minimises** — not `data − model`.



