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
integral equations. Lukas FEM is source-code reference material for clean-room
assembly patterns, not the API model.

When MATLAB classes are introduced, they should preserve the same feeling:
opening the class file should show the mesh, space, operator, or block tree in
plain mathematical terms. Hidden performance caches, clever vectorization, and
large opaque helper layers are secondary to student understanding.

## Basic API

```matlab
addpath("S:\MATLAB\Gypsilab\matlab_api");

m = volFemBem("mesh.vol");
uh = h1(m);       % H1 P1 tetrahedra
ah = hcurl(m);    % HCurl Nedelec0 tetrahedra
jh = rwg(m);      % boundary RWG0 dofs
```

## H-matrix Teaching Path

```matlab
H = educationalLaplaceHMatrix(m);
y = educationalHMatrixMatvec(H, ones(H.size(2), 1));
stats = educationalHMatrixStats(H);
```

The implementation is intentionally explicit:

- cluster tree
- admissibility test
- dense near-field blocks
- SVD low-rank far-field blocks
- recursive block-tree matvec

## Acoustic Teaching Path

```matlab
op = educationalAcousticSingleLayer(m, [], "Wavenumber", 10.0);
p = op.apply(ones(size(op.matrix, 2), 1));
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
K = lowFrequencyStableHelmholtzKernel(x, y, "Wavenumber", 1e-6);
K.singleLayerLaplace
K.singleLayerCorrection
K.singleLayer
```

## Optimization Teaching Path

```matlab
A = [1 0; 0 1; 1 1];
b = [1; 2; 3];
gate = educationalQuadraticLeastSquares(A, b, "Initial", zeros(2, 1));
gate.gradientCheck
```

This keeps the optimization lesson small and inspectable: the objective,
gradient, least-squares solution, initial point, finite-difference check, and
normal-equation residual are all returned in one struct. Optuna operation is
left to the official Optuna MCP server; this repository only keeps the MATLAB
teaching gate.

FEM/BEM trace fitting uses the same readable optimization style:

```matlab
m = volFemBem("mesh.vol");
fit = educationalFemBemTraceLeastSquares(m, boundaryValues, "Tikhonov", 1e-6);
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
m = volFemBemModel("mesh.vol");
manifest = educationalFemBemCouplingManifest(m);
manifest.checks
```

The manifest records the `.vol` mesh id, surface mesh id, trace artifact id,
volume space, surface space, formulation id, BEM kernel family, and
`traceRowIdentity`.  The trace operator also records the same row identity, so a
notebook can read which FEM node and BEM node each boundary row represents
without inferring it from sparse nonzeros.  A wrong expected kernel family is
`needs_attention` even when the trace matrix itself is one-hot and dimensionally
valid.

The same manifest records `boundaryRowIdentity`: each boundary triangle row is
bound to its triangle nodes, boundary number, boundary name, and adjacent
tetrahedron.  This keeps boundary-condition labels attached to the exact
surface triangle rows before BEM kernels, flux rows, or coupling notebooks reuse
the trace.

Touchstone rows used as MATLAB optimization objectives or constraints should
first carry port/option-line metadata, then pass a design-frequency sweep check
and the solver-ready row preflight:

```matlab
meta = educationalTouchstonePortMetadata(["P1" "P2"], ...
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
grid = educationalTouchstoneDesignFrequencyGrid([0.95 0.99 1.01 1.05], 1.0, ...
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
row = educationalTouchstoneSolverReadyPreflight(0.05, 0.8*exp(-1i*pi/18), ...
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
pat = educationalFarfieldPatternMetadata([0 90 180], [0 90], ...
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

This keeps frequency, RI/MA/DB format, reference impedance, passivity,
reciprocity, and S11 match quality visible before the row enters an
optimization notebook.

For the regularization trade-off, keep the path visible:

```matlab
path = educationalFemBemTikhonovPath(m, boundaryValues, [0; 1e-3; 1e-1; 1]);
path.checks
```

Larger Tikhonov weights should reduce the FEM unknown norm and increase the
boundary trace residual.  This is the readable L-curve-style lesson before
larger FEM/BEM inverse problems.

When a noise norm is known, the same path can be selected by Morozov's
discrepancy principle:

```matlab
choice = educationalFemBemMorozovDiscrepancy(m, boundaryValues, ...
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
box = educationalBoxConstrainedLeastSquares(A, b, [0; 0], [1; 3]);
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
eq = educationalTouchstoneEquivalentCircuit(0, 0, ...
    "S12", 0, "S22", 0, "Z0", 50, "ComparisonZ0", 75);
eq.pi.yShunt1   % 1/Z0 for the matched isolated teaching row
eq.t.zSeries1   % Z0 for the same row
```

This helper deliberately keeps `Z0` visible.  The same normalized S-parameter
row gives different equivalent admittance/impedance values when the Touchstone
reference impedance changes, so `data_format` and `Z0` must be recorded before
pi/T extraction.

## Mesh Policy

- Mesh handoff is Netgen `.vol`.
- Surface elements must be triangles.
- Volume elements must be tetrahedra.
- Only first-order tri/tet is supported in this lane.
- Quad, hex, wedge, pyramid, and hidden splitting are intentionally rejected.
- Third-party or high-order mesh evidence should be checked before `.vol`
  intake with a readable manifest:

```matlab
meshGate = educationalMeshImportQualityManifest(["tri" "tri"], ["tet"], ...
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

For acoustic FEM/BEM cases, COMSOL's acoustic FEM/BEM workflow is an important
internal secondary reference. It is used to understand coupling conventions,
radiation/open-boundary modeling, and low-frequency behavior. COMSOL references
must stay private and must use the already-running LiveLink MATLAB session.
