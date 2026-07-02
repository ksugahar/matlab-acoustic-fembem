# Class API Refactor (2026-07)

Status: executed 2026-07-03. This document records the design, the old-to-new
name map, and the verification protocol of the readable-class refactor.

## Goal

Raise `matlab_api` readability to the Gypsilab source level, following
`READABLE_CLASS_STYLE.md`:

- one mathematical object per class (mesh, space, trace, kernel, H-matrix,
  coupled system)
- short methods, explicit properties, no hidden caches
- assembly helpers live next to their class, Gypsilab `openMsh/mshXxx` style
- the `educational` prefix is dropped: the whole repository is educational,
  so the prefix carried no information
- teaching report gates stay plain struct-returning functions; they are
  reports, not mathematical operators, so a class would add ceremony only

## Directory layout

```
matlab_api/
  mesh/      readVolTriTet.m  VolMesh.m  SurfaceMesh.m
  fem/       H1Space.m  Nedelec0Space.m
  bem/       SurfaceP1Space.m  RwgSpace.m  TraceOperator.m
  kernel/    HelmholtzKernel.m
  hmatrix/   HMatrix.m
  acoustic/  AcousticSingleLayer.m  lowFrequencyHelmholtzReport.m
             helmholtzKernelManifest.m
  model/     FemBemModel.m  femBemCouplingManifest.m  (+ manifest helpers)
  gates/     (teaching report gates, prefix-free)
```

Consumers add the whole API with `addpath(genpath(<repo>/matlab_api))`.

## Name map

Core objects (function + struct -> classdef):

| Old | New |
|---|---|
| `readVolTriTet(file)` | `readVolTriTet(file)` (unchanged parser, moved to `mesh/`) |
| `volFemBemModel(file)` | `VolMesh(file)` + `SurfaceMesh(mesh)` + identity on the objects |
| `volFemBem(file)` | `FemBemModel(file)` |
| `addFirstOrderFemBemSpaces(model)` | folded into `FemBemModel` (space catalog property) |
| `buildFirstOrderTopology(model)` | folded: tet edges -> `Nedelec0Space`, tri edges -> `RwgSpace`, trace maps -> `TraceOperator` / `FemBemModel` |
| `h1(model)` | `H1Space(mesh)` |
| `hcurl(model)` | `Nedelec0Space(mesh)` |
| `rwg(model)` | `RwgSpace(surface)` |
| `assembleLukasP1Stiffness(model, c)` | `H1Space.stiffness(c)` |
| `assembleLukasP1Mass(model, c)` | `H1Space.mass(c)` |
| `assembleNedelec0TetMatrices(model, c)` | `Nedelec0Space.matrices(c)` (also `.mass(c)`, `.curlCurl(c)`) |
| `assembleGypsilabP1SurfaceMass(model)` | `SurfaceP1Space.mass()` |
| `assembleFirstOrderFemBemTrace(model)` | `FemBemModel.assembleTrace()` -> `TraceOperator` |
| `lowFrequencyStableHelmholtzKernel(x,y,...)` | `HelmholtzKernel(x,y,...)` |
| `educationalLaplaceHMatrix(...)` | `HMatrix(...)` |
| `educationalHMatrixMatvec(H,x)` | `HMatrix.matvec(x)` (and `H * x`) |
| `educationalHMatrixStats(H)` | `HMatrix.stats()` |
| `educationalAcousticSingleLayer(...)` | `AcousticSingleLayer(...)` |

Teaching report gates (renamed functions, same struct contracts, moved to
`gates/` unless noted):

| Old | New |
|---|---|
| `educationalFemBemCouplingManifest` | `femBemCouplingManifest` (`model/`, split into helpers) |
| `educationalFemBemNormalFluxSignReport` | `femBemNormalFluxSignReport` |
| `educationalFemBemTraceMassReport` | `femBemTraceMassReport` |
| `educationalFemBemTraceLeastSquares` | `femBemTraceLeastSquares` |
| `educationalFemBemTikhonovPath` | `femBemTikhonovPath` |
| `educationalFemBemLcurveCorner` | `femBemLcurveCorner` |
| `educationalFemBemMorozovDiscrepancy` | `femBemMorozovDiscrepancy` |
| `educationalQuadraticLeastSquares` | `quadraticLeastSquares` |
| `educationalBoxConstrainedLeastSquares` | `boxConstrainedLeastSquares` |
| `educationalGeometricIntegratorEnergyReport` | `geometricIntegratorEnergyReport` |
| `educationalMeshImportQualityManifest` | `meshImportQualityManifest` |
| `educationalMqCoulombGaugePostprocessPackage` | `mqCoulombGaugePostprocessPackage` |
| `educationalLowFrequencyHelmholtzKernelManifest` | `helmholtzKernelManifest` (`acoustic/`) |
| `lowFrequencyHelmholtzTeachingReport` | `lowFrequencyHelmholtzReport` (`acoustic/`) |
| `educationalTouchstonePortMetadata` | `touchstonePortMetadata` |
| `educationalTouchstonePortMatch` | `touchstonePortMatch` |
| `educationalTouchstoneEquivalentCircuit` | `touchstoneEquivalentCircuit` |
| `educationalTouchstoneDesignFrequencyGrid` | `touchstoneDesignFrequencyGrid` |
| `educationalTouchstoneSolverReadyPreflight` | `touchstoneSolverReadyPreflight` |
| `educationalFarfieldPatternMetadata` | `farfieldPatternMetadata` |
| `educationalFarfieldGainDirectivityEfficiency` | `farfieldGainDirectivityEfficiency` |
| `educationalFarfieldLobeNotebookHandoff` | `farfieldLobeNotebookHandoff` |

Error identifiers and `result.kind` strings follow the new function names
(`educational*:` ids become the new prefix-free ids, kind strings drop
`educational_`). Status strings that describe mathematical state
(`mesh_ready`, `vol_ready_first_order_h1_hcurl_rwg`,
`operators_ready_first_order_h1_hcurl_rwg_trace`) are unchanged.

Provenance-named struct containers disappear from the public surface:
`model.lukas.geo` becomes `VolMesh` data, `model.gypsilab.*` becomes
`SurfaceMesh` data. The optional real-Gypsilab interop is the explicit
`SurfaceMesh.gypsilabMsh()` method (errors when Gypsilab `msh` is not on the
path) instead of a silent constructor branch.

## Consumer surface

The 102 example scripts are 4-to-6-line wrappers around `runExampleById` /
`runMeshTopologyExample` and contain no direct API calls; they are untouched.
The real consumers rewritten to the class API:

- `tests/*.m` (11 files)
- `validation/verifyCatalogCase.m` and the mesh-topology helpers
- `examples/runExampleById.m` (unchanged unless it touches API names)

## Verification protocol

1. Baseline: full `run_tests.m` green before any edit (2026-07-03, EXIT=0).
2. Each stage adds the new classes next to the old functions and is gated by
   a throwaway parity probe in `C:\temp` (new assembly output must equal the
   old output exactly: same sparse triplets, same dense matrices, same ids).
3. Consumers switch to the new API, the old functions are deleted in the same
   change (single canonical API, no wrapper layer), and the full suite runs
   green again.
4. Final: full `run_tests.m` + 100-case batches green, validation log
   recorded under `S:\MATLAB\_crossval`.
