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

## Mesh Policy

- Mesh handoff is Netgen `.vol`.
- Surface elements must be triangles.
- Volume elements must be tetrahedra.
- Only first-order tri/tet is supported in this lane.
- Quad, hex, wedge, pyramid, and hidden splitting are intentionally rejected.

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
