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
