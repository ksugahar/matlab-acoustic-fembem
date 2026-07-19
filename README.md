# MATLAB Acoustic FEM/BEM

Readable MATLAB prototypes for acoustic FEM/BEM: tetrahedral volume FEM,
triangular boundary BEM, Laplace/Helmholtz kernels, FEM/BEM coupling,
fluid-structure interaction, and Lubich convolution-quadrature time-domain
solvers.

This repository is an educational solver laboratory.  It favors short,
inspectable MATLAB classes and reproducible validation artifacts over
production speed, and it is designed to be **driven through an MCP server**
rather than read as a hand-written API manual.

## How To Use This Repository (MCP-first)

**There is no generated API reference, and this README does not enumerate one.**
The public surface is deliberately just two things:

1. the short, readable MATLAB classes under `matlab_api/` -- open the file when
   you need the exact signature or convention; each class is meant to be read in
   one sitting; and
2. an **MCP server** that exposes the solver's gates, cross-validation, and
   teaching knowledge as tools.

The MCP layer is the intended entry point for AI agents, and the fastest way for
a human to drive a workflow without memorizing signatures.  It runs on the
official **MathWorks MATLAB MCP Server** as the execution substrate; this
repository supplies only the acoustic FEM/BEM domain tools, through the
`+acoustic_fembem` package and the extension file
[`mcp/extensions/acoustic-fembem-tools.json`](mcp/extensions/acoustic-fembem-tools.json):

| MCP tool | Purpose |
|----------|---------|
| `acoustic_fembem_knowledge` | Serve compact FEM/BEM/CQ teaching knowledge and conventions |
| `acoustic_fembem_vol_mesh_summary` | Summarize a `.vol` mesh (counts, bbox, orientation) with viewer guidance |
| `acoustic_fembem_acoustic_gate` | Run acoustic gates vs analytic references (`focus_adjoint`, `radiation_force`, `thrust_adjoint`) |
| `acoustic_fembem_crossval_gate` | `.vol`-backed cross-validation against NGSolve / `ngsolve.bem` references |
| `acoustic_fembem_check_result_manifest_file` | Validate a result manifest (versions, dates, schema/convention ids) |
| `acoustic_fembem_cq_time_grid` | Validate convolution-quadrature time and Laplace grids |
| `acoustic_fembem_cq_response_reality` | Check CQ conjugate symmetry and real time response |
| `acoustic_fembem_soft_sphere_cq_causality` | Gate the causal soft-sphere CQ response |
| `acoustic_fembem_adjoint_scaling` | Gate acoustic adjoint scaling behavior |
| `acoustic_fembem_hmatrix_scaling` | Gate readable ACA/H-matrix scaling behavior |

Setup and the full tool contract live in [mcp/README.md](mcp/README.md).  Install
the official MathWorks MATLAB MCP Server, then start it with the extension:

```powershell
--extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

For direct MATLAB use -- reading the classes, running a solve by hand -- see
[Quick Start](#quick-start) below, but treat the readable classes as the
reference, not a separate API document.

## What Is Included

- Netgen `.vol` intake and **pure-MATLAB `.vol` read/write** (`readVolTriTet`,
  `writeVol`, `icosphereVol`, `structuredBoxVol`) -- no PDE Toolbox required.
- Volume FEM spaces: H1 P1 tetrahedra and HCurl/Nedelec0 tetrahedra.
- Boundary BEM spaces: scalar P1 traces and RWG0 edge basis functions.
- FEM/BEM trace maps with explicit row identity.
- Dense Galerkin Laplace and Helmholtz BEM kernels, with low-frequency-stable
  Helmholtz splitting.
- **Curved (isoparametric) panel BEM** (`CurvedPanelQuadrature`,
  `curvedSingleLayerDirichletSolve`): Lagrange curve order 1/2/3 geometry with a
  P1 density -- the readable "curved P1" element that removes the O(h^2)
  straight-panel faceting error (~10-200x closer to the analytic sphere at the
  same mesh; a single `Projection` knob recovers the flat lane exactly). Curve
  order is the geometry lever, not fes order. An optional netgen path
  (`curvedQuadratureFromNetgen` + `tools/export_curved_boundary_nodes.py`)
  consumes a real curved mesh via a convention-free node companion, so the same
  accuracy gain works for general geometry.
- A readable educational H-matrix with **reference-guarded ACA+** low-rank
  far-field blocks -- rows and columns are sampled on the fly and the dense far
  block is never formed.
- Acoustic scattering, FEM/BEM coupling, FSI, radiation-force, and adjoint
  teaching gates.
- Lubich convolution-quadrature time-domain BEM (single-layer), including the
  drum-roll / scatterer demonstration.
- Reverse-mode adjoint gates for acoustic wavefront and radiation-force
  optimization.

The detailed solver story, including the validation ladder and acoustic/FSI
notes, lives in [docs/TEACHING_LADDER.md](docs/TEACHING_LADDER.md).

## Quick Start

```matlab
repoRoot = "path/to/matlab-acoustic-fembem";
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

## MATLAB `.vol` Read/Write (no PDE Toolbox)

Meshes can be created, written, and read back entirely in MATLAB:

```matlab
icosphereVol("ball.vol", Radius=1, Subdivisions=2, ...
    BoundaryName="outer", MaterialName="air");   % closed volume ball
mesh = VolMesh("ball.vol");                       % round-trips through readVolTriTet
```

`writeVol` emits the same first-order Netgen `mesh3d` contract used by the
FEM/BEM classes (surface triangles + optional volume tetrahedra, boundary and
material names).  For toolbox users, MATLAB PDE Toolbox export
(`writePdeMeshVol`) remains available as an optional path; it is not required to
build `.vol` meshes.

## Mesh Policy

The mesh handoff is Netgen `.vol`.

- Surface elements must be triangles.
- Volume elements must be tetrahedra.
- Only first-order tri/tet meshes are accepted.
- Quad, hex, wedge, pyramid, curved, and high-order records fail loudly.
- No hidden splitting or automatic topology conversion is performed.

This strict policy keeps the MATLAB teaching layer simple and makes trace
maps between volume and boundary spaces easy to inspect.

## Mesh Visualization

Netgen is the preferred interactive viewer for native `.vol` files.  MATLAB
keeps the lightweight figure/script and MCP preflight path:

```matlab
plotVolMesh(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_coarse.vol"));
summary = acoustic_fembem.vol_mesh_summary("unit_sphere_coarse.vol");
```

Use Netgen when you need to inspect the native mesh GUI, boundary labels, and
volume cells interactively.  Use `plotVolMesh` when a MATLAB script or
interactive session needs a quick boundary preview, and use the MCP tool
`acoustic_fembem_vol_mesh_summary` for headless counts, bounding boxes,
orientation, and viewer guidance.

## Repository Layout

```text
matlab_api/
  mesh/      .vol read/write (readVolTriTet/writeVol/icosphereVol/structuredBoxVol),
             PDE Toolbox exporters, preview plotter, VolMesh/SurfaceMesh
  fem/       H1Space, Nedelec0Space
  bem/       SurfaceP1Space, RwgSpace, TraceOperator
  kernel/    HelmholtzKernel
  hmatrix/   HMatrix (reference-guarded ACA+)
  acoustic/  acoustic reports, single-layer and convolution-quadrature helpers
  model/     FemBemModel and coupling manifests
  gates/     readable report gates for optimization, ports, far fields, traces

+acoustic_fembem/   MATLAB MCP gate package (domain tools over the MathWorks server)
mcp/                MCP extension file, verifier, and the tool contract
fixtures/           committed .vol meshes used by tests and validation
examples/           worked validation example entry points
validation/         cross-code / reference generators and runners
tests/              fast MATLAB unit-test lane
validation_test/    heavier numerical validation lane (analytic + cross-code)
docs/               design notes and the detailed teaching ladder
tools/              repository utilities
```

## Run Tests

Two lanes, matching the fast-regression / validation split:

```matlab
repoRoot = "path/to/matlab-acoustic-fembem";
run(fullfile(repoRoot, "run_tests.m"))             % fast lane   (tests/)
run(fullfile(repoRoot, "run_validation_test.m"))   % validation  (validation_test/)
```

`tests/` is the quick developer loop of fast unit checks.  `validation_test/`
holds the heavier numerical validation: analytic (Faran) FSI comparisons, the
convolution-quadrature drum roll, resonance sweeps, and cross-code references.

Focused validation examples can also be run directly:

```matlab
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

runMeshTopologyExample("GYP-001");
cases = validationCatalog();
```

## Requirements

- MATLAB R2021a or later.
- No required external solver for the fast unit-test lane.
- Optional: the official MathWorks MATLAB MCP Server (with the MATLAB Agentic
  Toolkit) to drive the repository through MCP tools.
- Optional: NGSolve / `ngsolve.bem` for regenerating cross-code reference
  artifacts.
- Optional: MATLAB PDE Toolbox, only for the `writePdeMeshVol` export path;
  `.vol` meshes can otherwise be built with `writeVol` / `icosphereVol`.

## Validation Status

- Fast lane (`tests/`) and validation lane (`validation_test/`) are both green.
- Acoustic and BEM checks against analytic references (e.g. Faran elastic-sphere
  scattering), committed NGSolve / `ngsolve.bem` reference artifacts, and MATLAB
  fixtures.
- Mesh/topology checks comparing MATLAB `.vol` intake with NGSolve `.vol`
  intake.

Generated logs are written under the platform temporary directory and are not
required for a basic checkout.

## MATLAB Execution Policy

This repository does not use MATLAB Live Editor documents as the default
interface.  Durable workflows should be ordinary `.m` functions/scripts, MCP
tools that print compact JSON, and result manifests with versions, dates,
timing, schema ids, and convention ids.

## Adjoint Optimization

This repository includes readable reverse-mode adjoint gates for optimization.
The key idea is to differentiate the equation, not the linear solver: one
transpose/adjoint solve gives the gradient of a scalar objective with respect
to many source amplitudes or design variables.

Current gates include:

- acoustic focusing: phased-array amplitudes maximize `|p(target)|^2`;
- acoustic radiation force: the force is a quadratic form of the linear field;
- elastic-bead thrust: FSI solve + Brillouin stress + Wirtinger gradient.

The MCP entry point `acoustic_fembem_acoustic_gate` exposes these as
`focus_adjoint`, `radiation_force`, and `thrust_adjoint`.

## Relationship To Gypsilab

Gypsilab is the style reference: compact notation, readable operators, and
MATLAB code close to the mathematics.  This repository does not vendor the
Gypsilab source tree.  Optional cross-checks against a separately installed
Gypsilab checkout are kept behind explicit validation helpers.

Because Gypsilab is GPL-3.0 software and this project intentionally follows
that ecosystem, this repository is released under GNU GPL v3.0.

## License

GNU GPL v3.0.  See [LICENSE](LICENSE).
