# Netgen .vol FEM/BEM coupling design for Lukas FEM + Gypsilab

Status: design plus first-order mesh/topology prototype, not a validated
coupled solver yet.

## Scope

This lane uses Cubit/Coreform Netgen `.vol` as the only mesh handoff format.
The parser accepts only:

- surface triangles
- volume tetrahedra

It rejects quad, hex, pyramid, and wedge records.  There is no hidden splitting
or tetrahedralization after export.

## Current working part

The current converter can build a shared-node and shared-edge coupling view:

- `geo`: Lukas-style volume FEM mesh fields
- `gypsilab`: Gypsilab `msh(vtx, elt, col)` surface mesh fields
- `trace`: one-based node ids shared by FEM volume and BEM boundary
- `topology.h1`: H1 P1 volume node ids
- `topology.hcurl`: HCurl/Nedelec0 volume edge ids and edge orientation signs
- `topology.scalarBem`: compacted boundary P1 node ids for Gypsilab-style BEM
- `topology.rwg`: RWG0 boundary edge dofs and triangle adjacency
- `topology.trace`: H1-to-scalar-BEM and RWG-to-HCurl trace maps

This is the mesh and first-order topology contract.  It does not yet assemble a
validated FEM/BEM linear system.

## First solver target

The API should already expose all first-order spaces needed for the intended
FEM/BEM lane:

- FEM: volume P1/H1 tetrahedra
- FEM: volume HCurl/Nedelec0 tetrahedra
- BEM: boundary P1 scalar traces
- BEM: boundary RWG0 triangle-pair currents/fluxes

The first numerical solve target is still scalar Laplace or magnetostatic
scalar-potential coupling because it isolates `.vol` trace maps and BEM sign
conventions.  The H(curl)/RWG topology must exist from the beginning so edge
orientation bugs are caught before any Maxwell kernel work.

The H(curl)/RWG case must not be treated as a simple node trace:

- Lukas FEM: Nedelec edge unknowns in the volume
- Gypsilab: RWG unknowns on the boundary surface
- coupling map: boundary oriented edges, not only boundary nodes

## Readable MATLAB API (class form, 2026-07)

Keep the user-facing path as short as Gypsilab:

```matlab
m = FemBemModel("model.vol");
m.h1          % H1Space, P1 tetrahedra
m.hcurl       % Nedelec0Space, HCurl tetrahedra
m.rwg         % RwgSpace, boundary RWG0 dofs
m.trace       % TraceOperator, H1 -> boundary P1
m = m.assemble();
```

One class per mathematical object (see `READABLE_CLASS_STYLE.md` and the name
map in `CLASS_API_REFACTOR.md`); the model stays boring and inspectable:

```matlab
m.mesh.vtx / m.mesh.tet / m.mesh.tri     % VolMesh, parsed .vol + identity
m.mesh.traceNodeIds
m.surface.vtx / m.surface.tri            % SurfaceMesh, compact boundary
m.surface.rowIdentity                    % boundary-condition row identity
m.spaceCatalog()
m.operators                              % after assemble()
```

Lukas FEM is source-code reference material for clean-room assembly style.  It
should not dictate the public MATLAB API.  The public API should remain
Gypsilab-like: small names, explicit spaces, and readable classes whose
properties expose the mathematical data directly.

## Readability over performance

This MATLAB lane is not meant to catch NGSolve or NGSolve.BEM on speed.  Its
job is to teach the solver architecture:

- how `.vol` tetrahedra become H1 and HCurl volume spaces
- how boundary triangles become scalar BEM and RWG spaces
- how traces connect FEM volume dofs to BEM boundary dofs
- how dense near-field BEM blocks and low-rank far-field blocks form an
  H-matrix
- how a block tree matvec represents what high-performance BEM libraries do
  with more engineering

The code should therefore prefer short functions, explicit structs, direct
linear algebra, and readable variable names.  Avoid clever vectorization if it
hides the mathematics.  Once the idea is correct and tested here, NGSolve,
NGSolve.BEM, or a compiled backend can carry the performance work.

If MATLAB classes are introduced, they should follow the same rule.  A class
should expose one mathematical object, such as a mesh, space, kernel, block, or
coupled system.  Its source should remain close enough to the formula that a
student can read it without first learning a private optimization framework.
See `READABLE_CLASS_STYLE.md`.

## Cross validation ladder

1. Mesh identity:
   - boundary triangles use only volume mesh node ids
   - Gypsilab `msh(vtx, tri, col)` uses compacted boundary ids with an explicit
     map back to one-based volume ids
   - boundary area and volume match radia-mcp `.vol` metrics

2. First-order topology identity:
   - H1 volume nodes match `.vol` point ids
   - HCurl volume edges and signs match radia-mcp `.vol` topology
   - RWG boundary manifold edges map back to HCurl boundary edges

3. Scalar manufactured solution:
   - unit tetra/sphere-like mesh
   - analytic harmonic potential
   - FEM volume residual and BEM boundary residual separately checked
   - FEM half landed 2026-07: `laplaceDirichletSolve` (interior Dirichlet
     partition-eliminate solve, P1 linear patch test locked in
     `tests/testLaplaceDirichletSolve.m`); the BEM-side residual check is
     still ahead

4. Exterior Laplace sphere:
   - capacitance or Dirichlet-to-Neumann map
   - compare Gypsilab dense BEM, Gypsilab hmx BEM, and radia-ngsolve

5. FEM/BEM coupled scalar open-boundary solve:
   - interior finite domain with exterior BEM boundary
   - compare boundary trace and integrated flux with radia-ngsolve reference

6. H(curl)/RWG Maxwell or magnetostatic vector coupling:
   - only after edge-orientation coupling and scalar signs are tested

## Gypsilab hmx performance expectation

Gypsilab has useful H-matrix pieces:

- `hmx`
- ACA compression
- QR-SVD recompression
- hierarchical subdivision
- H-matrix LU/Cholesky/LDL wrappers
- FEM/BEM EFIE/CFIE non-regression examples

But it should be treated as an algorithmic MATLAB prototype, not assumed to be
HACApK-class performance.  The implementation is MATLAB-recursive, and some
H-matrix algebra paths are incomplete.  HACApK remains the production-grade
compiled backend target for large problems.

For this repository, that is a feature rather than a flaw: Gypsilab's value is
that it makes H-matrix BEM unusually readable in MATLAB.  The local
`HMatrix` class intentionally copies that teaching
quality, not its full feature set.

## MEX candidates

MEX should wait until the scalar coupling is correct.  Then prioritize:

1. Near-field singular and regularized BEM triangle-pair quadrature.
2. Green-kernel block sampling used by ACA.
3. H-matrix matvec traversal and low-rank leaf application.
4. H-matrix factorization/preconditioner kernels if direct hmx solves are kept.
5. Lukas element assembly only if profiling shows it dominates.

Do not MEX the `.vol` parser, trace maps, or small sparse incidence builders
first.  They are readability-critical and unlikely to dominate.
