# Getting started

## Install

```bash
pip install pyrump
```

## Quick start

```bash
pyrump                         # the interactive shell, from any directory
```

```
Your wish? cd examples          /* cd to the examples directory from repo root  */
Your wish? get Sample1.rbs      /* read a spectrum and its metadata             */
Your wish? plot 1               /* plot the data from buffer 1                  */
Your wish? sim                  /* Get into simulation processor                */
SIM Command: layer 1            /* select layer 1                               */
SIM Command: comp Ta 1 N 1      /* set  composition of TaN                      */
SIM Command: thick 250 A        /* set thickness in A or nm                     */
SIM Command: comp Ta 1 N 1      /* sele       */
SIM Command: layer 2            /* sele       */
SIM Command: comp Si 1 N 1      /* sele       */
SIM Command: return
Your wish? compare              /* data vs simulation, with residuals */
```

Buffer 0 is always the simulation and recomputes itself when the sample or
the active buffer's parameters change — there is no "simulate" command,
exactly as in the original. See [Interactive shell](../manual/shell.md) for
the full session and command set.

Or drive it as one-off batch commands:

```bash
pyrump simulate sample.lcm --energy 2.0 --beam 4He -o out.rbs
pyrump fit sample.lcm measured.rbs --vary thickness:0 --window 190 226
pyrump plot measured.rbs --compare out.rbs -o comparison.png
pyrump convert measured.rbs measured.dat
```

See [CLI reference](../manual/cli.md) for every option, or [Python API](../dev/python-api.md)
to call the library directly. For a walked-through example of identifying and
fitting a real sample, see [Worked examples](examples.md).
