# Netgen .vol FEM/BEM coupling design for Lukas FEM + Gypsilab

Status: the cross-validation ladder below is climbed through stage 6
(2026-07): interior Dirichlet FEM solve, exterior Galerkin single-layer BEM
(sphere capacitance, cross-checked against the real Gypsilab AND
ngsolve.bem), the Johnson-Nedelec coupled FEM/BEM scalar open-boundary
solve (unit-ball source against the analytic radial solution), and the
H(curl)/RWG vector coupling layer (trace identity, RWG operators,
magnetostatic sphere gate). The acoustic lane opened 2026-07 with stage 7
(Helmholtz single-layer exterior solve, 3-way validated: analytic /
this repo / ngsolve.bem). A full vector transmission SOLVE (eddy-current
FEM/BEM with the vector Calderon operators) is intentionally beyond this
teaching lane.

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
   - landed 2026-07: Galerkin first-kind solve `singleLayerDirichletSolve`
     over `GalerkinSingleLayer` (test-side Gauss quadrature + all-pairs
     analytic Laplace panels `laplacePanelIntegrals` + quadrature for the
     smooth Helmholtz correction). Coarse unit sphere: C = 12.205 vs
     4*pi (geometry-faceting-dominated, -2.9%); same-mesh cross-check vs
     real Gypsilab integral+regularize: operator 1.1e-4, capacitance
     1.4e-5 (`validation/verifyGalerkinAgainstGypsilab.m`)
   - radia-ngsolve leg landed 2026-07: same-mesh, same-continuous-P1-space
     dense reference from NGSolve's `ngsolve.bem` (Sauter-Schwab,
     intorder 16, self-converged to 1e-8; exported by
     `validation/exportNgsolveBemReference.py`, checked by
     `validation/verifyGalerkinAgainstNgsolve.m` +
     `tests/testNgsolveBemCrossCheck.m` from the committed .mat).
     Measured: mass 1e-16, V 3.8e-4 / K 3.4e-3 at gss 7, capacitance
     1.6e-4 (coarse) / 4.3e-5 (fine). Conventions match exactly
     (outward-normal principal-value K: K[1] = -1/2 to 1e-9; interior
     H1 rows exactly zero). Honest finding: the gss-3 agreement with
     Gypsilab (1.1e-4) was same-test-quadrature error cancellation; the
     true assembly error vs the converged reference is 7e-3 at gss 3,
     4e-4 at gss 7 (V), matching the internal refinement study

5. FEM/BEM coupled scalar open-boundary solve:
   - interior finite domain with exterior BEM boundary
   - compare boundary trace and integrated flux with radia-ngsolve reference
   - landed 2026-07: `femBemCoupledSolve` (Johnson-Nedelec pair
     `A u - T' M lambda = F`, `(1/2 M - K) T u + V lambda = 0`) over
     `GalerkinDoubleLayer` (outward-normal principal value; sphere
     spectral gates K[1] = -1/2 exact to 6 digits, K[Y_1] = -1/6 to
     0.3%). Unit-ball source f = 1 against the analytic radial solution
     u = 1/2 - r^2/6: trace mean -2.3%, flux conservation
     int lambda dS = -int f dV to 1e-3, exterior potential 3.5%
     (all geometry-faceting dominated, improving coarse -> fine mesh).
     The kernel-sign conventions were locked numerically (Gauss check
     -4*pi inside, two spherical-harmonic BIE modes), not on paper
   - radia-ngsolve leg 2026-07: the coupled system's BEM ingredients
     (V and the principal-value K) are pinned operator-level against
     ngsolve.bem on the same meshes (see stage 4); the coupled solve
     itself keeps the analytic ball solution as its reference, which
     is stronger than a second numerical code on that geometry

6. H(curl)/RWG Maxwell or magnetostatic vector coupling:
   - only after edge-orientation coupling and scalar signs are tested
   - landed 2026-07 (the Maxwell-precursor layer):
     - the FEEC trace identity, pointwise to machine precision on every
       boundary edge: n x (tangential trace of the volume Nedelec0 edge
       function) = -signOut * triEdgeSigns * sigma_pm * f_RWG / l_e, and
       its dof-level form `RwgSpace.rotatedTraceMap` (RWG coefficients of
       n x u|_Gamma are C * alpha, exactly)
     - RWG basis machinery (`basisAtQuadrature`, exact 3-point `gram`)
     - `RwgSingleLayer`: the Galerkin vector single layer (static EFIE /
       partial-inductance kernel), assembled from the SAME analytic P1
       panel integrals (RWG components are affine per triangle - no new
       singular math), with `vectorPotentialAt`
     - magnetostatic gate: K = z_hat x n on the unit sphere (uniformly
       magnetized sphere) reproduces A = (1/3) z_hat x x inside and the
       dipole field outside (coarse 3.1%, fine 1.3%, faceting dominated)
     - found and fixed en route: readVolTriTet adjacentTetIndices had an
       ismember argument-order bug that silently returned tet 1 for every
       boundary triangle (orientation signs stayed correct for convex
       meshes; the adjacency values were wrong on multi-tet meshes)
   - intentionally NOT here: the full vector transmission solve
     (eddy-current FEM/BEM needs the vector double layer and
     hypersingular Calderon operators - production work for
     NGSolve.BEM, not this readable teaching lane)

7. Acoustic Helmholtz exterior solve (the acoustic-simulator lane opener):
   - landed 2026-07, THREE-way validated (analytic / this repo /
     ngsolve.bem) on the sphere fixtures:
     - `singleLayerDirichletSolve` gained the `Wavenumber` option: V_k by
       the existing GalerkinSingleLayer Helmholtz path, and a k-aware
       `potentialAt` with the SAME split (analytic Laplace panels + smooth
       expm1 correction), so the k -> 0 limit is exact (measured: density
       7e-10, probe potential 5e-11 against the Laplace solve)
     - analytic references in `matlab_api/acoustic`: `acousticPointSource`
       (interior source -> exterior reproduction, the EXACT gate) and
       `softSphereScattering` (partial-wave series; soft-BC residual on
       the true sphere 2.4e-12, truncation tail reported)
     - ngsolve.bem reference artifacts (committed .mat, k = 0.5 and 2.0,
       intorder 16 self-converged to ~3e-9, NODAL boundary data, plus
       NGSolve's own GetPotential probe values via auxiliary probe meshes):
       `exportNgsolveBemHelmholtzReference.py` /
       `verifyHelmholtzAgainstNgsolve.m` / `testHelmholtzScattering.m`
     - measured: operator V reldiff 1.2e-3..8.7e-3 (vs conjugate 0.68..1.3
       => e^{+ikr} pinned both sides); probe cross-code 2.7e-4..6.2e-3
       while both codes sit 1-10% from the true sphere (faceting,
       improving ~x2.5 coarse -> fine) - the two codes agree 10-30x
       tighter than either matches the analytic sphere, so the analytic
       deviation is geometry, not implementation
     - taught caveat: the first-kind V_k equation is singular at interior
       Dirichlet eigenvalues (unit sphere: kR = pi); k = 0.5 / 2.0 sit
       safely away. CHIEF / Burton-Miller (and the Helmholtz double layer
       K_k = K_0 + smooth correction) are the next acoustic rungs, then
       rigid scattering and the interior Helmholtz FEM + coupled
       transmission problem.

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
