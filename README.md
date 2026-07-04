# CAE-AI MATLAB FEM/BEM

Readable MATLAB prototypes for learning the FEM/BEM ideas behind Gypsilab,
NGSolve, and NGSolve.BEM.

This repository is an educational solver laboratory.  It favors short,
inspectable MATLAB classes and reproducible validation artifacts over
production speed.  The main target is a student-friendly FEM/BEM codebase that
shows how tetrahedral volume FEM, triangular boundary BEM, trace maps, acoustic
kernels, and small optimization gates fit together.

## What Is Included

- Netgen `.vol` mesh intake for first-order triangle/tetrahedron meshes.
- Volume FEM spaces: H1 P1 tetrahedra and HCurl/Nedelec0 tetrahedra.
- Boundary BEM spaces: scalar P1 traces and RWG0 edge basis functions.
- FEM/BEM trace maps with explicit row identity.
- Dense Galerkin Laplace and Helmholtz BEM kernels.
- Low-frequency-stable Helmholtz kernel splitting.
- A readable educational H-matrix matvec.
- Acoustic scattering, FEM/BEM coupling, FSI, radiation-force, and adjoint
  teaching gates.
- 100 validation examples organized under `examples`.

The detailed solver story, including the validation ladder and acoustic/FSI
notes, lives in [docs/TEACHING_LADDER.md](docs/TEACHING_LADDER.md).

## Relationship To Gypsilab

Gypsilab is the style reference: compact notation, readable operators, and
MATLAB code close to the mathematics.  This repository does not vendor the
Gypsilab source tree.  Optional cross-checks against a separately installed
Gypsilab checkout are kept behind explicit validation helpers.

Because Gypsilab is GPL-3.0 software and this project intentionally follows
that ecosystem, this repository is released under GNU GPL v3.0.

## Requirements

- MATLAB R2021a or later.
- No required external solver for the fast unit-test lane.
- Optional: NGSolve / NGSolve.BEM for regenerating cross-code reference
  artifacts.
- Optional: a separately installed Gypsilab checkout for the explicit
  Gypsilab comparison helper.

## Quick Start

```matlab
repoRoot = "path/to/caeai-matlab-fembem";
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "validation"));

m = FemBemModel(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol"));

m.mesh          % VolMesh: nodes, tetrahedra, boundary triangles
m.surface       % SurfaceMesh: compact boundary view
m.h1            % H1Space: volume P1
m.hcurl         % Nedelec0Space: volume edge space
m.scalarBem     % SurfaceP1Space: boundary P1
m.rwg           % RwgSpace: boundary RWG0
m.trace         % TraceOperator: FEM boundary trace
```

As a first numerical solve:

```matlab
surface = VolMesh(fullfile(repoRoot, "fixtures", "mesh_topology", ...
    "unit_sphere_coarse.vol")).boundary();

sol = singleLayerDirichletSolve(surface, ones(size(surface.vtx, 1), 1));
sol.totalCharge
sol.potentialAt([3 0 0])
```

## Mesh Policy

The mesh handoff is Netgen `.vol`.

- Surface elements must be triangles.
- Volume elements must be tetrahedra.
- Only first-order tri/tet meshes are accepted.
- Quad, hex, wedge, pyramid, curved, and high-order records fail loudly.
- No hidden splitting or automatic topology conversion is performed.

This strict policy keeps the MATLAB teaching layer simple and makes trace
maps between volume and boundary spaces easy to inspect.

## Repository Layout

```text
matlab_api/
  mesh/      .vol parser, VolMesh, SurfaceMesh
  fem/       H1Space, Nedelec0Space
  bem/       SurfaceP1Space, RwgSpace, TraceOperator
  kernel/    HelmholtzKernel
  hmatrix/   HMatrix
  acoustic/  acoustic reports and single-layer helpers
  model/     FemBemModel and coupling manifests
  gates/     readable report gates for optimization, ports, far fields, traces

fixtures/    committed .vol meshes used by tests and examples
examples/    100 validation example entry points
validation/  cross-code/reference generators and validation runners
tests/       fast MATLAB unit tests
docs/        design notes and the detailed teaching ladder
```

## Run Tests

From MATLAB:

```matlab
repoRoot = "path/to/caeai-matlab-fembem";
run(fullfile(repoRoot, "run_tests.m"))
```

Focused validation examples can also be run directly:

```matlab
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

runMeshTopologyExample("GYP-001");
cases = validationCatalog();
```

## Validation Status

The public validation campaign currently has:

- `100 / 100` example scripts present.
- `100 / 100` catalog entries marked verified.
- Mesh/topology checks comparing MATLAB `.vol` intake with NGSolve `.vol`
  intake.
- Acoustic and BEM checks against analytic references, committed NGSolve /
  NGSolve.BEM reference artifacts, and MATLAB fixtures.

Generated logs are written under a local sibling `_crossval` directory and are
not required for a basic checkout.

## MATLAB MCP Integration

The companion repository
[caeai-matlab-mcp](https://github.com/ksugahar/caeai-matlab-mcp) exposes a
MATLAB MCP domain layer around selected gates in this FEM/BEM repository.
The MCP repository treats the official MathWorks MATLAB MCP Server as an
external runtime and keeps domain behavior in versioned MATLAB functions.

## License

GNU GPL v3.0.  See [LICENSE](LICENSE).
