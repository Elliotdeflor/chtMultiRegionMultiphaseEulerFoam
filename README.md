# chtMultiRegionMultiphaseEulerFoam

Conjugate heat transfer (CHT) solver for Euler-Euler flows with any number of
phases, for ESI OpenFOAM v2606. Each fluid region is a `multiphaseSystem`
solved with the equations of `reactingMultiphaseEulerFoam`; solid regions,
region handling and time-step control follow `chtMultiRegionTwoPhaseEulerFoam`.
Solid regions are optional, so the same solver also runs fluid-only cases with
a heat-flux or fixed-temperature wall.

builds and runs on v2606.


---

## 1. Install

### Requirements

- ESI OpenFOAM **v2606** (other versions untested: the solver borrows source
  files from the installed `chtMultiRegionFoam` and `reactingMultiphaseEulerFoam`).
- The compiler and `wmake` that come with that OpenFOAM installation.

### Build

```bash
source /usr/lib/openfoam/openfoam2606/etc/bashrc        # adjust to your install
cd <path>/chtMultiRegionMultiphaseEulerFoam              # this folder
wmake 2>&1 | tee log.wmake
which chtMultiRegionMultiphaseEulerFoam                  # -> $FOAM_USER_APPBIN/...
ls -l $(which chtMultiRegionMultiphaseEulerFoam)         # check the timestamp is new
```

The folder can live anywhere (the include paths in `Make/options` are
absolute). The executable always goes to `$FOAM_USER_APPBIN`.

### Things that break the build

| Symptom | Cause | Fix |
| --- | --- | --- |
| `'rhoReactionThermo' was not declared` from `chtMultiRegionFoam/fluid/createFluidFields.H` | `fluid/createFields.H` is missing, so the stock single-phase file is used instead (see "Include order" below) | Restore `fluid/createFields.H` (two lines, see the file table) and run `wclean && wmake` |
| `'ls' was not declared` or similar nonsense at line 1 of a file | Something was typed into a source file in the editor | Remove it |
| Build seems to succeed but runs behave as before | Two copies of the source; the other one was built last | Keep one source folder; check `ls -l $(which chtMultiRegionMultiphaseEulerFoam)` |

After adding or removing files, run `wclean` before `wmake` so the dependency
lists are rebuilt.

### Include order (why one file must keep a stock name)

`#include "X.H"` resolves to the **first** `X.H` found: the including file's
own folder, then `.`, then the `-I` folders in `Make/options` order. Nothing
warns about other files with the same name. Local `./fluid` and `./solid` come
first, so local files override upstream ones.

All local files have unique names except `fluid/createFields.H`, which must
keep its stock name because OpenFOAM's `postProcess.H` includes
`createFields.H` by name. If it is missing, the build silently picks up
`chtMultiRegionFoam/createFields.H` and fails much later with unrelated errors.

---

## 2. What is in the folder

| File | Purpose |
| --- | --- |
| `chtMultiRegionMultiphaseEulerFoam.C` | Main loop: time step, outer correctors over all regions, optional `energyCoupling` loop, write |
| `fluid/createFields.H` | **Required override** (stock name): `#include "createMultiphaseFluidFields.H"` and `#include "createSolidFields.H"` |
| `fluid/createMultiphaseFluidFields.H` | One `multiphaseSystem`, `p_rgh`, `gh`, `ghf`, `pMin`, `pimpleControl` per fluid region |
| `fluid/setRegionMultiphaseFluidFields.H` | Binds region `i` to the names the stock multiphase equation files use |
| `fluid/readMultiphaseFluidPIMPLEControls.H` | Per-region `faceMomentum`, `partialElimination`, `nEnergyCorrectors`, `frozenFlow`, `turbOnFinalIterOnly` |
| `fluid/multiphaseMultiRegionCourantNo.H` | Courant number: max over mixture flux, each phase flux and each phase pair's relative flux |
| `fluid/solveMultiphaseFluid.H` | One fluid outer iteration, in stock order; equation files (`UEqns.H`, `EEqns.H`, `pEqn.H`, `YEqns.H`) come from `reactingMultiphaseEulerFoam` |
| `fluid/readOuterPIMPLEControls.H`, `fluid/storeOldMultiphaseFluidFields.H` | Outer-loop controls and previous-iteration storage (`p_rgh` only) |
| `solid/*` | Solid energy, diffusion number: copies of the two-phase CHT solid files |
| `derivedFvPatchFields/turbulentTemperatureMultiphaseRadCoupledMixed` | N-phase fluid-solid coupling condition (not yet verified in a run) |
| `derivedFvPatchFields/turbulentTemperatureTwoPhaseRadCoupledMixed` | Unchanged copy of the two-phase condition, so two-phase CHT cases run as they are |

---

## 3. Setting up a case

### 3.1 Layout

Multi-region layout, as for any CHT solver:

```
case/
  constant/
    regionProperties          regions ( fluid (water) solid (heater) );
    g                         shared by all regions (top level, not in a region)
    water/  heater/           per-region physics and polyMesh
  system/
    controlDict               application chtMultiRegionMultiphaseEulerFoam;
    fvSolution                top-level PIMPLE: nOuterCorrectors (see 3.5)
    fvSchemes
    water/  heater/           per-region fvSchemes, fvSolution, setFieldsDict, ...
  0.orig/
    water/  heater/
  Tsat_*.csv, other tables    if phaseProperties refers to them by relative path
```

**CHT off:** `regions ( fluid (water) solid () );` and an ordinary wall
condition (heat flux or temperature) instead of the coupled one.

### 3.2 Converting a two-phase case

Starting from a `chtMultiRegionTwoPhaseEulerFoam` or `reactingTwoPhaseEulerFoam`
case, every item below was needed for `solidQuenching2D`:

1. `system/controlDict`: `application chtMultiRegionMultiphaseEulerFoam;`
2. `constant/<fluid>/phaseProperties`: two-phase system type -> multiphase
   equivalent (e.g. `thermalPhaseChangeTwoPhaseSystem` ->
   `thermalPhaseChangeMultiphaseSystem`). To list valid types, set `type xxx;`
   and run once.
3. `phaseProperties`: add `interfaceCompression ();` (empty = none). It is
   mandatory in the multiphase system.
4. Turbulence: `mixtureKEpsilon` is two-phase only. Use per-phase `kEpsilon`
   or `kOmegaSST`, and validate against a two-phase run with the **same**
   per-phase model.
5. **Every phase fraction needs real boundary conditions.** The two-phase
   solver derives the second phase as 1 - alpha and ignores its BCs, so
   tutorials often leave them as `calculated`. The multiphase solver
   transports every phase. At outlets use, for each phase,
   ```
   type        inletOutlet;
   phi         phi.<phase>;
   inletValue  uniform <value that may flow back in>;
   value       uniform <same>;
   ```
   In `solidQuenching2D`, `alpha.liquid` at the top boundary was `calculated`
   with value 1, which let liquid flow in through the outlet.
6. `setFieldsDict`: set **every** phase fraction, including a default for each.
7. Volume-fraction numerics (see 3.4).

### 3.3 Temperature conditions on a coupled wall

- Only one phase carries the coupled condition, normally the continuous
  liquid. Fluid side, e.g. `0.orig/water/T.liquid`:
  ```
  <patch>
  {
      type        compressible::turbulentTemperatureMultiphaseRadCoupledMixed;
      region      fluid;
      Tnbr        T;
      kappaMethod phaseSystem;
      // phases   (liquid gas);     // optional; default = all phases
      value       uniform 300;
  }
  ```
  Solid side, `0.orig/heater/T`:
  ```
  <patch>
  {
      type        compressible::turbulentTemperatureMultiphaseRadCoupledMixed;
      region      solid;
      Tnbr        T.liquid;         // the field that carries the condition
      kappaMethod solidThermo;
      value       uniform 300;
  }
  ```
  `otherPhase` is not used (a warning is printed if present). Wall temperature
  = (K_s T_s + sum_k K_k T_k) / (K_s + sum_k K_k), with K_k = alpha_k kappaEff_k / delta.
- **Every other phase's temperature** on that patch copies the wall value:
  ```
  <patch> { type copiedFixedValue; sourceFieldName T.liquid; value uniform 300; }
  ```
- Two-phase cases may keep `compressible::turbulentTemperatureTwoPhaseRadCoupledMixed`
  (with `otherPhase`); it is compiled into this solver.
- Heat-rate printout for conservation checks, in `system/controlDict`:
  ```
  DebugSwitches { compressible::turbulentTemperatureMultiphaseRadCoupledMixed 1; }
  ```

How boiling reaches the wall: `alphatWallBoilingWallFunction` back-calculates
each phase's wall `alphat` so that kappaEff (T_w - T_cell)/delta gives that
phase's share of the boiling heat flux. A **negative** liquid `alphat` at the
wall is normal in film boiling; it also makes the wall flux very sensitive to
small changes near the wall.

### 3.4 Volume fractions and time step

The multiphase system solves every phase fraction separately and rescales them
to sum to 1 at the end of each step. With settings inherited from two-phase
cases the per-step drift can be large and bias the wall heat flux: in
`solidQuenching2D` it gave +35% wall heat flux until fixed with:

```
// system/<fluid>/fvSolution
solvers
{
    "alpha.*"
    {
        nAlphaCorr      3;
        nAlphaSubCycles 3;      // 5 made the gas energy equation diverge
        nLimiterIter    10;
    }
    "(e|h|k|epsilon).*" { ...; maxIter 100; }
}
// system/controlDict
maxCo   0.25;
```

The ESI wall-boiling tutorials take the other route: default alpha settings
with `maxCo 0.05`. Either way, **monitor the phase sum** (section 4).

### 3.5 PIMPLE controls: where each one is read

| Entry | Read from | Notes |
| --- | --- | --- |
| `nOuterCorrectors` | **top-level** `system/fvSolution` | The region's own value is ignored. Copy it to the top level when converting a single-region case. |
| `nCorrectors`, `nNonOrthogonalCorrectors` | `system/<fluid>/fvSolution` | |
| `nEnergyCorrectors`, `faceMomentum`, `partialElimination`, `frozenFlow` | `system/<fluid>/fvSolution` | |
| `turbOnFinalIterOnly` | `system/<fluid>/fvSolution` | Default `true`, as in the stock solver: turbulence corrected on the final outer iteration only |

With `nOuterCorrectors` > 1, the `energyCoupling` loop re-solves fluid energy
(flow frozen) and solid energy on the first outer iteration. Near CHF,
where the heat flux changes steeply with wall temperature, use several outer
correctors and check that the interface heat flux converges within each step.

### 3.6 Other case files

- **Function objects** in `controlDict` need `region <name>;`, and
  `postProcess` needs `-region <name>`. Without it they fail with
  "failed lookup of region0" or silently look at the wrong mesh.
- `constant/g` and data tables referenced by relative path (e.g.
  `Tsat_water_1_2bar.csv`) sit at the case top level.
- Field relaxation: only `p_rgh` keeps a previous-iteration copy. This matters
  only if a case field-relaxes another field with `nOuterCorrectors` > 1.

---

## 4. Running and monitoring

```bash
./Allrun.pre                                  # mesh, regions, setFields
nohup ./Allrun > Allrun.log 2>&1 &            # survives closing the terminal
```

- `runApplication` **skips any step whose `log.*` already exists**. To rerun:
  stop the solver (`pgrep -af EulerFoam`), then `rm -f log.*`, remove
  `processor*` and old time folders.
- Tutorials whose `Allrun` edits `controlDict` (start time, end time) keep a
  `controlDict.orig`; restore it before rerunning.
- Parallel: `decomposePar -allRegions`, run, `reconstructPar -allRegions`.
- Logs contain non-text bytes: use `grep -a`.
- Keep case paths free of spaces; OpenFOAM refuses a working directory with a
  space in it.

**Checks every run:**

```bash
L=log.chtMultiRegionMultiphaseEulerFoam
grep -a -c "FOAM FATAL\|nan" $L                                  # 0
grep -a "Phase-sum" $L | awk '{n++; if($8<0.99)a++; if($8<0.9)b++} \
    END{print n" checks; min<0.99: "a"; min<0.9: "b}'             # few or none below 0.9 after start-up
```

Instability signature: gas energy residual near 1 for many iterations and the
temperature limiter clipping many cells, starting in cells with almost no gas.

---

## 5. Differences from the stock solvers

| Behaviour | Here | Stock `reactingMultiphaseEulerFoam` |
| --- | --- | --- |
| Outer loop | Plain loop over all regions (multi-region design) | `pimple.loop()` |
| Early exit on converged outer residuals | No | Optional (`outerCorrectorResidualControl`) |
| First time step | Not adjusted when Courant and diffusion numbers are zero (flow at rest, no solids) | Adjusted before the first step |
| Courant number | Also includes relative phase fluxes | Mixture and phase fluxes |
| Turbulence correction | Final outer iteration only (`turbOnFinalIterOnly`) | Same |

Because of the first two time-step differences, adaptive runs do not match
the stock solver step for step. For code-equivalence tests use a fixed time
step (`adjustTimeStep no`).

The stock `chtMultiRegionTwoPhaseEulerFoam` corrects turbulence on **every**
outer iteration, so it differs from `reactingTwoPhaseEulerFoam` whenever
`nOuterCorrectors` > 1.

---\

## 6. Known limits

- `alphatWallBoilingWallFunction` pairs one liquid with one vapour (`otherPhase`)
  and its vapour side assumes liquid fraction = 1 - vapour fraction. Fine for
  two phases; with more phases at the wall it is approximate until replaced by
  the wall-boiling library wrapper.
- Energy coupling between regions is segregated (no single fluid-solid matrix);
  convergence relies on outer correctors.
- No local time stepping.
- If a run stops with an unknown turbulence model, the multiphase turbulence
  library is not loaded: find it with `ls $FOAM_LIBBIN | grep -i turbulence`
  and add it to `libs ( ... );` in `system/controlDict`.

Licence: GPL v3 or later (derived from OpenFOAM).
