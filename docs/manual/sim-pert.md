### SIM and PERT

`SIM` edits the sample description; `PERT` fits it. Both are sub-levels with
their own prompt, and — as in the original — a command the sub-level does not
recognise is passed out to the RUMP level, which returns you there automatically.
`SIM <command>` and `PERT <command>` also work as one-shots from the top level
-- e.g. `PERT GET usual.pert GO` loads a saved fit setup and runs it in a
single line, handy in a macro.

```
Your wish? sim
SIM Command: thick 500 A
SIM Command: composition Si 1 /
SIM Command: next
SIM Command: thick 2000 A
SIM Command: composition Au 1 /
SIM Command: show
SIM Command: save mysample.lcm
SIM Command: return
```

SIM's sample-definition commands are the *same code* that parses `.lcm` files
(`SampleEditor` in `pyrump/script/lcm.py`), so what you type and what the file
holds cannot drift apart. `SIM SAVE` writes RUMP's own format.

#### SIM commands

| Command | Effect |
| --- | --- |
| `?` / `HELP` | list the SIM commands |
| `RETURN` / `ABORT` / `QUIT` | return to the RUMP level |
| `LAYER <n>` | move to layer *n* |
| `NEXT` | move to the next layer, opening one if needed |
| `OPEN` | insert a blank layer above the current one |
| `DELETE` / `CLOSE` | remove the current layer |
| `RESET` | reset the sample to empty space |
| `SHOW` | display the sample description |
| `STATUS` | summarize layers, maxpth, straggle, multiple |
| `THICKNESS <v> <unit>` | this layer's thickness |
| `COMPOSITION El n [El n …] /` | this layer's stoichiometry |
| `SPECIES El n [El n …] /` | the impurity species an `EQUATION` blends toward |
| `EQUATION <name> <params…>` | depth-profile equation for this layer |
| `EQLIST` | list the known equation names |
| `FUZZ <amount> <steps>` | roughen the interface above this layer |
| `SUBLAYER <n>` | force a sublayer count |
| `STHICKNESS <v> <unit>` | or set the thickness of each sublayer |
| `MAXPTH [<v>]` | default sublayer thickness, 10¹⁵ at/cm² -- shows the current value with no argument |
| `STRAGGLE <v>` | Bohr straggling multiplier |
| `ABSORBER <n>` | first *n* layers are a dead layer/window, not sample |
| `MULTIPLE <v>` | multiple-scattering amount |
| `GET <file>` | read a sample description from a `.lcm` file |
| `SAVE <file>` | write the sample description to a `.lcm` file |
| `DENSITY [pattern]` | list known thickness units, or matching `density.tab` compounds |
| `SPLOT` | overlay the simulation on the current plot |
| `COMPARE` | active buffer vs. the simulation, with residuals |

```
SIM Command: layer 1
SIM Command: thick 500 A
SIM Command: composition Si 1 /
SIM Command: density ito       /* look up how "ITO" resolves as a Thick unit */
```

PERT selects what may vary and over which channels, then `GO`:

```
Your wish? pert
PERT Command: window 355 375      /* compare only here            */
PERT Command: norm 140 200        /* rescale data to remove dose error */
PERT Command: thickness 1         /* vary layer 1's thickness     */
PERT Command: go
```

```
  Fitting SC0338: Si [5000/cm2] - Au [200/cm2]

  fit took 0.34 s
  reduced chi-square 1.2849 on 20 dof
  12 evaluations, Both `ftol` and `xtol` termination conditions are satisfied.
  data scaled by 0.99441 over the norm window
  layer 1 thickness                    299  +/- 0.3801   (was 200)

  SC0338: Si [5000/cm2] - Au [299/cm2]
```

`GO` now opens with "Fitting `<ID>`: `<structure>`" *before* the fit, and
closes with the same `<ID>` (bare, never the buffer's raw path/label) and the
structure *after* -- `<ID>` is the active data buffer's own file stem (see
`REPORT`), so it's stable even if a WRASCII macro stamped a full path into the
buffer's name via `FILENAME`.

Fitted values are written back into the sample description, so `SIM SHOW` and
`SIM SAVE` reflect them. Two differences from the original: the data may be in
any buffer, not just buffer 1, and `MULTI` is the default because the solver is
a simultaneous least-squares fit (`SINGLE` loops one parameter at a time).

#### Reading the fit report

`reduced chi-square 1.2849 on 20 dof` is the standard Pearson chi-square per
degree of freedom for Poisson-counting data. `dof` is the number of channels
inside the active error `WINDOW`(s) minus the number of varying parameters,
floored at 1 (`chi_square()`, `fit/objective.py`); the reduced value is just
`chi-square / dof`. As a rule of thumb: close to 1 means the model explains
the data about as well as counting statistics allow; **much greater than 1**
points at a systematic mismatch the solver cannot fit away — a missing
layer, a bad calibration or geometry value, or a window that includes
non-Rutherford scattering (not modeled here — see the
[known-limitations note](../dev/about.md#milestones)) — rather than noise,
since noise alone averages out over hundreds of channels; **much less than
1** usually means the error window is too narrow or otherwise not exercising
real Poisson statistics. A residual plot with structure (a slope, a bump the
model tracks smoothly through) rather than a flat scatter around zero is
normally where the mismatch actually is.

`N evaluations, <status>` reports how many full spectrum simulations `GO`
ran, and whether the fit `converged` -- one of `scipy.optimize.least_squares`'s
tolerances (`xtol` or `ftol`, both set to RUMP's `EpsCrit`) was satisfied and
the solver settled on a minimum -- or `did not converge`, meaning it hit its
evaluation cap first. Either way the reported parameters are whatever the
solver last tried, so `did not converge` is a signal to check the starting
values or narrow the window rather than trust the numbers as a settled fit.

#### PERT commands

| Command | Effect |
| --- | --- |
| `?` / `HELP` | list the PERT commands |
| `RETURN` / `QUIT` | return to the RUMP level |
| `GO` | run the fit |
| `PARMS` | display the current selection and windows |
| `SHOW` | display the sample description (same as `SIM SHOW`) |
| `GET <file>` / `GET <file> GO` | replay a saved selection from a `.pert` file, replacing the current one, and optionally run the fit right after |
| `SAVE <file>` | save the current selection (windows + varying parameters) to a `.pert` file |
| `CLEAR` / `CLEAR <n>` | forget everything, or just the *n*th varying parameter (1-based, as numbered by `PARMS`) |
| `WINDOW lo hi` / `WINDOW clear` / `WINDOW clear <n>` | add / clear all / clear the *n*th error window, in channels (up to 10) |
| `NORMALIZE lo hi` / `NORMALIZE clear` | set / clear the normalisation window |
| `SINGLE` / `MULTI` | fit one parameter at a time / all together (default) |
| `VOLUME [off]` | print a line per model evaluation during `GO` |
| `AUTOCMP [off]` | run `COMPARE` automatically at the end of `GO` (default off) |
| `REPORT [off]` | after every `GO`, save `<sample>.report` (fit results, appended), `<sample>.pert` (`PERT SAVE`), `<sample>.lcm` (`SIM SAVE`) and `<sample>.png` (`FIGSAVE`) -- `<sample>` is the active data buffer's own file stem, so switching samples with `XEQ`/`GET` routes later fits to different files automatically (default off) |
| `THICKNESS <layer> [<min> <max>]` | vary a layer's thickness, optionally bounded |
| `COMPOSITION <layer> <El> [<min> <max>]` | vary one element's composition in a layer (must already be declared there) |
| `SPECIES <layer> <El> [<min> <max>]` | vary the `EQUATION` species composition (must already be declared there) |
| `EQUATION <layer> <n> [<min> <max>]` | vary equation parameter *n* |
| `MEV` / `FWHM` / `THETA` / `CORRECTION` / `STRAGGLE` `[<min> <max>]` | vary that beam, detector or sample parameter |
| `SLOPE` / `KEV/CH` `[<min> <max>]` | vary the calibration slope, keV per channel |
| `OFFSET` / `KEV(0)` `[<min> <max>]` | vary the calibration energy offset alone (e.g. a sample-charging shift) |
| `FUZZ` | not implemented — raises an error |
| `COMPARE` `[new]` | active buffer vs. the simulation, with residuals |

Every one of the above takes an optional trailing `<min> <max>` search
bound -- a real constraint the solver enforces, not just a label. Omit it and
the parameter keeps its default range; give it and it replaces that default
outright. `PARMS` echoes any bound in force, and `SAVE`/`GET` round-trip it.

> **Never vary every element's composition in a layer — always leave one
> fixed.** A layer's composition is stored as raw, unnormalized stoichiometric
> coefficients (e.g. `Co 1 Mn 0.333197 Si 0.324647`), and it's normalized by
> its own sum before it ever reaches the physics — both in the slab fill
> (`fractions = row / total`, `sim/slabs.py`) and in the Bragg mixing rule for
> matrix density (`ρ = (Σxᵢ/ρᵢ / Σxᵢ)⁻¹`, `atomic/density.py`), which cancels
> any common scale exactly. So scaling *every* composition in a layer by the
> same factor changes nothing in the simulated spectrum, for any factor — it's
> an exactly flat direction in the fit, not just a poorly-conditioned one. For
> an *N*-element layer, vary at most *N*-1 of its `COMPOSITION` selections and
> leave one as the fixed reference (conventionally the majority element, at
> its nominal coefficient); the others then fit as genuine ratios to it.
> Violating this shows up as a singular Jacobian and unusable uncertainties
> (`fit/lm.py`'s `_covariance` catches exactly this and reports `None`), or
> the solver parking on an arbitrary point wherever the bounds happen to stop
> it.

```
PERT Command: window 355 375   /* compare only here                  */
PERT Command: thickness 1      /* vary layer 1's thickness           */
PERT Command: mev              /* also vary the beam energy          */
PERT Command: go
```

