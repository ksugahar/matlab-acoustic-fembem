# MATLAB Gypsilab Education Layer

Readable MATLAB prototypes for learning FEM/BEM ideas behind Gypsilab,
NGSolve, and NGSolve.BEM.

This repository is a local, non-public lab repository. It is not a production
solver and it is not meant to compete with NGSolve on speed. Its job is to make
the mathematics readable:

- Cubit/Coreform Netgen `.vol` intake
- first-order H1 P1 tetrahedral FEM
- first-order HCurl/Nedelec0 tetrahedral FEM
- scalar boundary P1 BEM traces
- RWG0 boundary edge dofs
- FEM/BEM trace maps
- dense and low-rank BEM blocks
- educational H-matrix matvecs
- acoustic Helmholtz BEM kernels
- low-frequency-stable acoustic BEM kernels
- readable quadratic least-squares optimization gates
- readable box-constrained projected-gradient optimization gates

Gypsilab is the style reference: short notation, clear source, and enough
operator structure that a reader can connect MATLAB code to the boundary
integral equations. The classes follow `docs/READABLE_CLASS_STYLE.md`: one
mathematical object per class, public properties that expose the mathematics,
short methods, no hidden performance caches. The name map and layout of the
2026-07 class refactor are recorded in `docs/CLASS_API_REFACTOR.md`.

## Basic API

```matlab
addpath(genpath("S:\MATLAB\Gypsilab\matlab_api"));

m = FemBemModel("mesh.vol");
m.mesh          % VolMesh: vtx / tet / tri + labels + source identity
m.surface       % SurfaceMesh: compact boundary + row identity
m.h1            % H1Space (volume P1)
m.hcurl         % Nedelec0Space (volume edge)
m.scalarBem     % SurfaceP1Space (boundary P1)
m.rwg           % RwgSpace (boundary triangle-pair dofs)
m.trace         % TraceOperator: g = m.trace * u

[K, femDetail] = m.h1.stiffness();        % int grad(u).grad(v) dx
[M, bemDetail] = m.scalarBem.mass();      % int u v dS
[Mn, Cn, edgeDetail] = m.hcurl.matrices();

m = m.assemble();   % operators struct for the coupling manifest
```

Each class file is meant to be read: opening `H1Space.m` shows the P1
barycentric-gradient assembly loop, opening `TraceOperator.m` shows the
one-hot injection and its artifact identity, opening `RwgSpace.m` shows the
surface edge extraction and the RWG-to-HCurl oriented-edge map.

## Boundary Value Problem Teaching Path

The first BVP rung of the cross-validation ladder is the interior Laplace
Dirichlet solve (partition and eliminate, no BEM kernel yet):

```matlab
m = FemBemModel("mesh.vol");
g = ...;                             % one value per boundary trace node
sol = laplaceDirichletSolve(m, g);   % u_B = g,  K_II u_I = -K_IB g
sol.u
sol.interiorResidualNorm
sol.checks
```

The test suite locks the P1 patch test: for constant coefficient and linear
boundary data the discrete solution reproduces the linear potential exactly
(1e-12).

The exterior-BEM rung (ladder stage 4) is the Galerkin single-layer solve:

```matlab
mesh = VolMesh("fixtures/mesh_topology/unit_sphere_coarse.vol");
surface = mesh.boundary();
sol = singleLayerDirichletSolve(surface, ones(size(surface.vtx, 1), 1));
sol.totalCharge          % capacitance for g = 1: 12.205 vs 4*pi = 12.566
u = sol.potentialAt([3 0 0]);
```

`GalerkinSingleLayer` follows the real Gypsilab `integral + regularize`
split: the test integral uses Gauss quadrature (`SurfaceQuadrature`, 1/3/7
points), the singular Laplace kernel is integrated analytically over every
source triangle (`laplacePanelIntegrals`, Wilton-style closed forms verified
to machine precision against subdivision and polar references), and the
smooth low-frequency-stable Helmholtz correction goes through plain
quadrature. Cross-checked against the real Gypsilab on the same sphere mesh:
operator relative difference 1.1e-4, capacitance relative difference 1.4e-5
(`validation/verifyGalerkinAgainstGypsilab.m`). The remaining ladder rung is
the coupled FEM/BEM scalar open-boundary solve (stage 5).

## H-matrix Teaching Path

```matlab
H = HMatrix(m);                    % points, SurfaceMesh, or FemBemModel input
y = H * ones(H.shape(2), 1);       % or H.matvec(x)
stats = H.stats();
```

The implementation is intentionally explicit:

- cluster tree
- admissibility test
- dense near-field blocks
- SVD low-rank far-field blocks
- recursive block-tree matvec

## Acoustic Teaching Path

```matlab
op = AcousticSingleLayer(m, [], "Wavenumber", 10.0);
p = op.apply(ones(op.shape(2), 1));    % or op * q
```

This starts with the dense Helmholtz single-layer operator
`exp(1i*k*r)/(4*pi*r)` so the acoustic BEM equation stays visible. Internally
it uses the low-frequency-stable split

```matlab
G_k = G_0 + (G_k - G_0)
```

where the correction is evaluated with `expm1` or a Taylor series. This keeps
the `k -> 0` limit connected to the Laplace BEM kernel and mirrors the
low-frequency-stability discipline used by serious BEM codes. Once the kernel
and signs are understood, the same ideas can be compressed with an H-matrix or
moved to NGSolve.BEM.

For direct kernel inspection:

```matlab
K = HelmholtzKernel(x, y, "Wavenumber", 1e-6);
K.singleLayerLaplace
K.singleLayerCorrection
K.singleLayer
```

## Optimization Teaching Path

```matlab
A = [1 0; 0 1; 1 1];
b = [1; 2; 3];
gate = quadraticLeastSquares(A, b, "Initial", zeros(2, 1));
gate.gradientCheck
```

This keeps the optimization lesson small and inspectable: the objective,
gradient, least-squares solution, initial point, finite-difference check, and
normal-equation residual are all returned in one struct. Optuna operation is
left to the official Optuna MCP server; this repository only keeps the MATLAB
teaching gate.

FEM/BEM trace fitting uses the same readable optimization style:

```matlab
m = FemBemModel("mesh.vol");
fit = femBemTraceLeastSquares(m, boundaryValues, "Tikhonov", 1e-6);
fit.gradientCheck
```

The objective is explicit:
`0.5*(T*u-g)'*M*(T*u-g) + 0.5*alpha*(u'*u)`, where `T` is the
FEM-to-BEM trace and `M` is the boundary P1 mass matrix.  This is the smallest
bridge from `.vol` trace maps to optimization before moving to real acoustic
or electromagnetic FEM/BEM coupling.

Before coupled values are compared, keep the trace and operator identity as one
package:

```matlab
m = FemBemModel("mesh.vol");
manifest = femBemCouplingManifest(m);
manifest.checks
```

`femBemCouplingManifest` RECORDS the package (mesh id, surface mesh id, trace
artifact ids, volume space, surface space, formulation id, BEM kernel family,
`traceRowIdentity`); `femBemManifestChecks` VALIDATES the recorded report and
is readable as the list of contract lines the package must satisfy.  The trace
operator carries its own row identity, so a notebook can read which FEM node
and BEM node each boundary row represents without inferring it from sparse
nonzeros.  A wrong expected kernel family is `needs_attention` even when the
trace matrix itself is one-hot and dimensionally valid.

Boundary-condition identity is single-sourced on the `SurfaceMesh`: each
boundary triangle row is bound to its triangle nodes, boundary number,
boundary name, and adjacent tetrahedron (`m.surface.rowIdentity`), and the
manifest records it before BEM kernels, flux rows, or coupling notebooks reuse
the trace.

Touchstone rows used as MATLAB optimization objectives or constraints should
first carry port/option-line metadata, then pass a design-frequency sweep check
and the solver-ready row preflight:

```matlab
meta = touchstonePortMetadata(["P1" "P2"], ...
    "PortOrder", ["P1" "P2"], ...
    "NetworkParameter", "S", ...
    "DataFormat", "MA", ...
    "FrequencyUnit", "GHz", ...
    "Z0", 50);
meta.checks
```

The row values are not interpreted until port names/order, network parameter,
RI/MA/DB format, frequency unit, and `Z0` are explicit.  This keeps swapped
ports and missing reference impedance out of later optimization examples.

```matlab
grid = touchstoneDesignFrequencyGrid([0.95 0.99 1.01 1.05], 1.0, ...
    "FrequencyUnit", "GHz", ...
    "MaxRelativeSpacing", 0.03);
grid.lowerIndex
grid.upperIndex
grid.bracketGapRel
```

The nearest row is not enough when interpolation, equivalent-circuit fitting,
or group-delay fitting will reuse the sweep.  Keep the lower/upper bracket rows
and `bracketGapRel = bracket_gap / design_frequency` visible before accepting
the row-level evidence.

```matlab
row = touchstoneSolverReadyPreflight(0.05, 0.8*exp(-1i*pi/18), ...
    "S12", 0.8*exp(-1i*pi/18), ...
    "S22", 0.05, ...
    "Frequency", 1.0, ...
    "FrequencyUnit", "GHz", ...
    "DataFormat", "MA", ...
    "ReturnLossMinDb", 20, ...
    "VswrMax", 1.2);
row.checks
```

Far-field pattern rows used in antenna, EMC, or ngsolve.bem teaching notebooks
follow the same metadata-first rule:

```matlab
pat = farfieldPatternMetadata([0 90 180], [0 90], ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RequiredPhiValuesDeg", [0 90], ...
    "RowCount", 6);
pat.checks
```

This keeps frequency, theta/phi cuts, polarization basis, quantity unit, row
count, and power normalization visible before gain/directivity/RCS values are
used in optimization or FEM/BEM comparison.

For the regularization trade-off, keep the path visible:

```matlab
path = femBemTikhonovPath(m, boundaryValues, [0; 1e-3; 1e-1; 1]);
path.checks
```

Larger Tikhonov weights should reduce the FEM unknown norm and increase the
boundary trace residual.  This is the readable L-curve-style lesson before
larger FEM/BEM inverse problems.

When a noise norm is known, the same path can be selected by Morozov's
discrepancy principle:

```matlab
choice = femBemMorozovDiscrepancy(m, boundaryValues, ...
    [0; 1e-3; 1e-2; 1e-1; 1], noiseNorm);
choice.selectedAlpha
choice.checks
```

The selected row is the one whose weighted trace residual is closest to the
noise norm.  Keeping this rule visible helps students separate "choose the
smooth-looking L-curve corner" from "choose the alpha justified by the data
noise estimate."

Bound-constrained design variables use the same inspectable style:

```matlab
A = eye(2);
b = [2; -1];
box = boxConstrainedLeastSquares(A, b, [0; 0], [1; 3]);
box.activeLower
box.activeUpper
box.maxKktResidual
```

The projected-gradient step is intentionally visible.  At a lower bound the
gradient should be nonnegative, at an upper bound it should be nonpositive, and
free variables should have zero gradient.  This is the small constrained gate
used before larger FEM/BEM design optimization notebooks.

Touchstone / port-mode rows can also be inspected as readable linear algebra
before moving them into a circuit or BEM notebook:

```matlab
eq = touchstoneEquivalentCircuit(0, 0, ...
    "S12", 0, "S22", 0, "Z0", 50, "ComparisonZ0", 75);
eq.pi.yShunt1   % 1/Z0 for the matched isolated teaching row
eq.t.zSeries1   % Z0 for the same row
```

This helper deliberately keeps `Z0` visible.  The same normalized S-parameter
row gives different equivalent admittance/impedance values when the Touchstone
reference impedance changes, so `data_format` and `Z0` must be recorded before
pi/T extraction.

## API Layout

```
matlab_api/
  mesh/      readVolTriTet  VolMesh  SurfaceMesh  bemCollocationPoints
  fem/       H1Space  Nedelec0Space
  bem/       SurfaceP1Space  RwgSpace  TraceOperator
  kernel/    HelmholtzKernel
  hmatrix/   HMatrix
  acoustic/  AcousticSingleLayer  lowFrequencyHelmholtzReport
             helmholtzKernelManifest
  model/     FemBemModel  femBemCouplingManifest  femBemManifestChecks
  gates/     teaching report gates (touchstone*, farfield*, femBem*,
             quadraticLeastSquares, boxConstrainedLeastSquares, ...)
```

Classes are the mathematical objects; gates are plain struct-returning report
functions (they are lessons, not operators).

## Mesh Policy

- Mesh handoff is Netgen `.vol`.
- Surface elements must be triangles.
- Volume elements must be tetrahedra.
- Only first-order tri/tet is supported in this lane.
- Quad, hex, wedge, pyramid, and hidden splitting are intentionally rejected.
- Third-party or high-order mesh evidence should be checked before `.vol`
  intake with a readable manifest:

```matlab
meshGate = meshImportQualityManifest(["tri" "tri"], ["tet"], ...
    "Order", 1, ...
    "MinScaledJacobianBefore", 0.36, ...
    "MinScaledJacobianAfter", 0.75, ...
    "NegativeJacobianCountBefore", 2, ...
    "NegativeJacobianCountAfter", 0, ...
    "CadConnectivityRecorded", true, ...
    "CadComplianceRecorded", true, ...
    "BoundaryConformityTolerance", 1e-6, ...
    "MaxBoundaryDistance", 4e-7);
meshGate.checks
```

  This is the teaching bridge from high-order curvilinear mesh literature:
  CAD connectivity/conformity/compliance and scaled-Jacobian quality are
  recorded, while MATLAB still rejects high-order, quad, hex, prism, wedge, or
  pyramid topology instead of silently converting it.

## Run Tests

```matlab
run("S:\MATLAB\Gypsilab\run_tests.m")
```

Every test file also runs standalone (each carries a `setupOnce` that adds
`matlab_api` recursively to the path).

## 100-Case Validation Campaign

The `examples` directory is prepared for 100 radia-ngsolve cross-validation
cases, organized as 10 categories with 10 cases each. The catalog lives in:

```matlab
cases = validationCatalog();
```

A case can become `verified` only after:

- the MATLAB example script exists
- the radia-ngsolve reference exists
- the comparison tolerance is declared
- the run passes
- the validation log is recorded under `S:\MATLAB\_crossval`

The test suite checks that the catalog has exactly 100 unique cases, that all
verified cases have example scripts, and that validation logs are recorded.

Current progress:

- `100 / 100` example scripts present
- `100 / 100` verified
- verified groups: all 10 categories
- validation log:
  `S:\MATLAB\_crossval\gypsilab_mesh_topology_10of100_20260624.md`
  `S:\MATLAB\_crossval\gypsilab_remaining_90of100_20260624.md`
  `S:\MATLAB\_crossval\gypsilab_class_api_refactor_20260703.md`

For acoustic FEM/BEM cases, COMSOL's acoustic FEM/BEM workflow is an important
internal secondary reference. It is used to understand coupling conventions,
radiation/open-boundary modeling, and low-frequency behavior. COMSOL references
must stay private and must use the already-running LiveLink MATLAB session.
