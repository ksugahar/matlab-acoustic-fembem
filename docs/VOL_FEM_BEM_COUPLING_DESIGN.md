# Netgen .vol FEM/BEM coupling design for Lukas FEM + Gypsilab

Status: the cross-validation ladder below is climbed through stage 6
(2026-07): interior Dirichlet FEM solve, exterior Galerkin single-layer BEM
(sphere capacitance, cross-checked against the real Gypsilab AND
ngsolve.bem), the Johnson-Nedelec coupled FEM/BEM scalar open-boundary
solve (unit-ball source against the analytic radial solution), and the
H(curl)/RWG vector coupling layer (trace identity, RWG operators,
magnetostatic sphere gate). The acoustic lane opened 2026-07: stage 7
(Helmholtz single-layer exterior solve, 3-way validated: analytic /
this repo / ngsolve.bem), stage 8 (sonic-crystal chain, 4-leg incl. an
NGSolve volume-FEM ABC leg; free-space chain shows NO Bragg gap), and
stage 10 (the coupled acoustic FEM/BEM transmission solve, pure MATLAB,
3 gated cases). A full vector transmission SOLVE (eddy-current FEM/BEM
with the vector Calderon operators) is intentionally beyond this teaching
lane.

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
  student can read it without first learning a specialized optimization framework.
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

8. Sonic-crystal chain: multi-body scattering, 4-leg validated:
   - landed 2026-07: five sound-soft spheres (R = 0.3, d = 1.5,
     `soft_sphere_chain_5.vol` + committed generator). No new BEM
     machinery - the all-pairs assembly takes any number of closed
     surfaces; the fixed adjacentTetIndices/orientation path is what
     multi-body correctness rides on.
   - four legs at k = 1.0 and 2.0944 (= pi/d):
     exact interior point source inside the MIDDLE sphere (1.8-2.0e-2,
     gss 3); Foldy monopole multiple scattering (new
     `foldyPointScattering`; 2.7e-2 at k = 1.0, honest degradation to
     ~15% at kd = pi measured against the exact series); ngsolve.bem
     (operator 7.7-9.2e-3 at gss 3, probe cross-code 4e-4..1.1e-3,
     conjugate 0.59-0.87 => convention re-pinned); NGSolve volume FEM
     (order 2, 69k tets, Dirichlet spheres, first-order Sommerfeld ABC
     at R_out = 6: agrees 4.4-7.5e-2 - the methods-diverse leg).
   - PHYSICS FINDING (negative result, locked as a test): the sparse
     free-space chain has NO Bragg stop band at k d = pi. Insertion loss
     is a flat broadband sub-wavelength-attenuation plateau
     (3.47..3.94 dB over k = 0.6..3.0, spread < 0.5 dB; BEM and Foldy
     agree, FEM confirms at the two locked k). Soft spheres scatter with
     scattering length ~ R down to k -> 0, so the chain attenuates
     broadband instead of developing a lattice gap; free-space leakage
     suppresses the Bragg mechanism a duct would confine.

9. Sonic-crystal band gap, landed 2026-07 (NGSolve FEM reference,
   committed .mat, locked by tests/testDuctBandGap.m):
   - duct unit cell a = 1.0, d = 1.5, sound-soft sphere R = 0.3 (the
     stage-8 chain family), OCC z-periodic identification + Floquet
     phase exp(-i q d), order-2 complex H1, ArnoldiSolver
   - GATES: empty-lattice bands match the analytic duct dispersion
     (incl. transverse families) to 2.0e-4; the empty duct with one-way
     port Robin conditions is transparent to ~1e-4
   - RESULT: band 1 [2.31, 2.52], Bragg gap [2.52, 3.65]; finite 5-cell
     transmission shows stop (T ~ 1e-4 below band 1 - soft inclusions
     cut off the long-wavelength plane mode) / pass (T up to 0.27) /
     stop (T ~ 1e-3 in the gap), aligned with the Bloch bands
     (contrast > 100). Confinement creates what stage 8 proved free
     space cannot.
   - bug class caught en route (recorded): OCC face classification by
     AREA thresholds silently turned the duct into a Dirichlet cavity -
     classify geometrically (face center on a wall plane), and the
     empty-lattice gate is what caught it.

12. Adjoint automatic differentiation (wavefront synthesis), landed
   2026-07: reverse-mode AD through the rigid-scattering BEM solve.
   `acousticFocusAdjoint` designs phased-array complex amplitudes to focus
   `J = |u(target)|^2` behind a rigid scatterer; the field is affine in the
   amplitudes (`u = w p`, `w = S0 + lambda' M S` from one transpose solve
   `A' lambda = d0'`), so the whole gradient is ONE adjoint solve
   independent of the source count. Measured
   (`tests/testAcousticFocusAdjoint.m`): forward affine residual 4e-18,
   adjoint vs central FD 1.7e-10, gradient ascent focuses monotonically.
   The reusable `doubleLayerPotentialMatrix` (the trace -> field row) was
   extracted from rigidScatteringSolve so the solve, CHIEF, and the adjoint
   share one implementation. Sign trap recorded: the steepest-ascent
   direction of the real objective is the Wirtinger derivative
   `dJ/dconj(p) = 2 u conj(w)`, not `2 conj(u) w`; a near-zero start hides
   the flip (any step raises |u|^2), a nonzero start exposes it, and the
   adjoint-vs-FD gate on the real/imag gradients catches it regardless.
   The intensity objective becomes the net radiation force (stage 13); the
   same adjoint then yields dF/d(phases) for wavefront-synthesised thrust.

13. Acoustic radiation force (ultrasonic thrust / manipulation), landed
   2026-07: acousticRadiationForce computes the net time-averaged force -
   the SECOND-order functional of the linear field - as the control-sphere
   integral of the Brillouin radiation-stress tensor
   T_ij = [|p|^2/(4 rho c^2) - rho|v|^2/4] delta_ij + (rho/2) Re[v_i v_j*],
   v = grad p/(i omega rho), F = -oint T.n dS. The same post-processor
   takes the analytic series OR the BEM total field.
   - GATED formula-free: F is INDEPENDENT of the control radius to ~1e-10
     (T is divergence-free in the source-free fluid - the primary gate, no
     external formula), F_z > 0 (downstream), axisymmetric, Y_p(kR=2) =
     0.752 (King 1934); BEM force matches the series to ~5% (faceting).
   - Reynolds-stress SIGN was the bug the radius-independence gate caught:
     the divergence-free tensor has +(rho/2)Re[v_i v_j*], not minus (a
     minus gave control-radius spread of 400% incl. a sign flip; the
     axisymmetry gate passed either way, only radius-independence exposed
     it - lock the tensor by a formula-free consistency gate, not by a
     remembered sign).
   - PHYSICAL air @ 40 kHz (lambda 8.58 mm; kR = 2 <-> Rs = 2.73 mm):
     F_z = Y_p pi Rs^2 <E>, <E> = p_rms^2/(rho c^2): 120 dB -> 0.05 uN,
     140 dB -> 5 uN, 160 dB -> 0.5 mN - the acoustic-levitation /
     micro-manipulation regime (mm-bead weight ~10 uN near 150 dB).
   Adjoint dF/d(phases) - feed this force as the objective to the
   acousticFocusAdjoint machinery (the force is a quadratic form in the
   surface trace t, so the same one-adjoint-solve reverse mode applies)
   for wavefront-synthesised thrust optimization.

14a. Elastic sphere scattering - the FSI (fluid-structure) analytic
   reference, landed 2026-07: elasticSphereScattering is the Faran (1951)
   solid elastic sphere (interior compressional + shear potentials; 3x3
   boundary system - radial displacement continuity, radial stress =
   -pressure, zero shear stress). This is the reference for the acoustic
   fluid-structure INTERACTION case, where FEM/BEM coupling truly earns
   its keep: an elastic scatterer has INTERNAL resonances (compressional +
   shear) that rigid/soft spheres cannot show.
   - VALIDATED against two INDEPENDENT references: ShearSpeed -> 0
     reproduces the Anderson fluid sphere to 1.3e-10 (exact fluid limit),
     a stiff/heavy solid reproduces the rigid sphere (3.7e-3), and the 3x3
     with shear converges to the 2x2 fluid limit (1.3e-3).
   - SIGN BUG found + fixed: the radial-stress boundary term is
     -lambda kL^2 A j_l(kL R), NOT + (sigma_rr = -p). The flipped sign
     PASSED the stiff->rigid limit but failed mu->0 by 20-70% - a
     self-consistent limit (vs my own rigid reimpl) is NOT enough; only
     the INDEPENDENT Anderson reference exposed it. Gate on both.
   - Elastic radiation force (acousticRadiationForce on the elastic field):
     a lucite-like sphere (rho 1.15, cL 1.6, cT 0.9 vs fluid) has Y_p(kR)
     rising to 2.83 at kR = 3 (an internal resonance well above the rigid
     ~0.75 plateau) then falling, control-radius independent to 1e-10.
   Its coupled solve landed as stage 14b (below).

14b. Acoustic FSI COUPLED SOLVE, landed 2026-07: fsiCoupledSolve - the
   genuine acoustic FEM/BEM coupling. Vector P1 elasticity FEM interior
   (elasticityMatrices: per-tet 12x12 B'DB stiffness + consistent mass)
   coupled to the acoustic BEM exterior (V_k, K_k) through the FSI
   interface (dynamic: sigma(u).n = -p n; kinematic:
   (1/(rho_f w^2)) dp/dn = u.n). Coupled block system (the scalar
   Johnson-Nedelec femBemCoupledSolve with a VECTOR interior; unknowns
   u, p_s, q_s):
     [ K - w^2 M      G'            0   ] [u  ]   [-G' p_inc]
     [ -rho_f w^2 G   0             Mb  ] [p_s] = [-Mb q_inc]
     [ 0              1/2 Mb - K_k  V_k ] [q_s]   [ 0       ]
   with G = int_Gamma mu (n.v) the interface coupling, exterior
   representation p_s = D_k[p_s|Gamma] - S_k[q_s]. Also promoted the
   reusable singleLayerPotentialMatrix (density -> field, refactored out
   of singleLayerDirichletSolve for DRY).
   - VALIDATED against the analytic elastic sphere (stage 14a): the STIFF
     limit reproduces rigidSphereScattering to 5e-4 / 4.2e-3 / 2.2e-4 -
     the formulation gate (the interface + BEM coupling is exact,
     independent of interior resolution); the elastic field CONVERGES to
     elasticSphereScattering under refinement (coarse 25/55/26% -> fine
     7.5/16.7/7.7% at kR = 2). The remaining ~10% is the P1 interior
     elasticity (shear wavelength ~ sphere size at kR = 2 - low-order FEM
     dispersion), NOT the coupling.
   - grad/typo bug fixed (elasticityMatrices: grad = inv(C)(2:4,:), the
     last 3 rows of the P1 shape-coefficient inverse).
   This is the acoustic analogue of the lab's core EM coupling
   (nonlinear iron + open boundary; eddy currents + open boundary), where
   the BEM / Sommerfeld / spherical-DtN / high-order impedance boundary path
   plays the acoustic open-boundary role.  Do not call the acoustic boundary
   a Kelvin boundary.
   NEXT: higher-order interior (P2 elasticity) for the ~10% shear-
   resolution gap, and the adjoint dF/d(phases) with the ELASTIC (FSI)
   scatterer for elastic-bead thrust design.

11. Rigid scattering + irregular frequencies + CHIEF, landed 2026-07:
   - total-field equation (1/2 M - K_k) t = M g_inc
     (rigidScatteringSolve; mode-verified against the locked
     conventions: (1/2 - K_l) t_l = a_l j_l with t_l = a_l (i/x^2)/h_l')
   - ngsolve.bem K_k added to the sphere reference artifacts:
     operator 1.8-3.1e-3 (conjugate 0.12-1.4), ngbem-side rigid solve
     cross-code trace 4e-6..1e-4 (second-kind conditioning visible);
     both codes 3.8-13% from the rigid-sphere series (shared faceting)
   - irregular-frequency lesson LOCKED: at kR = pi the solution is 96%
     wrong while cond(A) ~ 29 stays benign (the faceted eigenvalue
     shifts) - only the analytic gate sees it; CHIEF interior
     null-field rows restore 8.0e-2. Burton-Miller (hypersingular)
     remains the declared production alternative, out of the teaching
     lane.

10. Coupled acoustic FEM/BEM transmission (pure MATLAB), landed 2026-07:
    - the Helmholtz double layer closed the operator set:
      K_k = K_0 (analytic panels) + smooth base*(exp(z)(1-z)-1)
      correction via HelmholtzKernel SourceNormals (which already carried
      the split); sphere gates K_k[Y_l] = 1/2 + 1i k^2 j_l(k) h_l'(k) to
      3-5e-3 (k=0.5) / 2-3e-2 (k=2), k->0 = Laplace K to 5e-28
    - femBemCoupledSolve gained Wavenumber / InteriorWavenumber /
      DensityRatio / IncidentAmplitude: FEM row
      (1/rho1)(A - k1^2 Mv) u - T'M sigma = F + T'G_inc, BIE row
      (1/2 M - K_k)Tu + V_k sigma = (1/2 M - K_k) g_inc, exterior
      representation -S_k + D_k. NO absorbing boundary anywhere - the BEM
      row is the exact radiation condition (the stage-8 FEM leg's ABC
      cost 4-7%; here the same physics needs none)
    - three gated cases (fine sphere -> finer unit_ball_maxh018, so the
      P1 (k h)^2 resolution limit is a locked CONVERGENCE assertion):
      (1) k->0 regression vs the verified Laplace coupled solve (1e-9);
      (2) acoustic invisibility k1=k0, rho1=rho0 - the exact null gate
      (interior == p_inc, u_s == 0): 4.1e-2/2.6e-2 -> 1.3e-2/7.8e-3;
      (3) Anderson fluid sphere c1/c0=0.7, rho1/rho0=1.2 vs the
      partial-wave series: interior 13% -> 4.4%, probes 22% -> 7.3%
    - fluidSphereScattering lesson: the naive per-mode 2x2 transmission
      match is catastrophically ill-conditioned at high l (RCOND 1e-23,
      h_l huge vs j_l tiny) - eliminate through the interior
      log-derivative so the invisible case gives A_l = 0 EXACTLY; and the
      series term count must key on the FARTHEST evaluation point
      (k0 * r_max), not the sphere radius

10b. Frequency-sweep time-domain bridge, landed 2026-07:
    - `volFemBemIfftResponse` reads the same tri/tet `.vol` mesh, assembles
      H1/P1 volume FEM and boundary P1 Galerkin BEM, solves
      `femBemCoupledSolve` over positive Helmholtz frequency bins, multiplies
      by a real pulse spectrum, mirrors the spectrum with Hermitian symmetry,
      and uses `ifft` to produce a real pressure time history at exterior
      observation points
    - this is not a periodic sine animation and not yet
      convolution-quadrature TD-BEM; it is the readable frequency-domain
      FEM/BEM + inverse-FFT route toward full P1 volume-FEM / P1 surface-BEM
      transient acoustics
    - first gate: `tests/testVolFemBemIfftResponse.m` on `unit_tetra.vol`
      verifies tri/tet intake, P1 FEM/P1 BEM identity, frequency solve
      residuals, Hermitian/iFFT realness, and nonzero pulse response

10c. Convolution-quadrature TD-BEM bridge, landed 2026-07:
    - `volTdBemConvolutionQuadrature` is the first Lubich CQ time-domain
      BEM rung. It uses the `.vol` boundary triangles as a P1 Galerkin BEM
      surface, samples the BDF generating function `delta(zeta)`, evaluates
      Laplace-domain retarded single-layer matrices `V(s)` at
      `s = delta(zeta)/dt` with `Re(s)>0`, solves `V(s) qhat = ghat`, and
      FFT-recovers the boundary density `q(t)` and exterior pressure
    - this is not a Helmholtz frequency-sweep/iFFT response; the CQ Laplace
      points come from the time-stepping formula. Current scope is exterior
      Dirichlet single-layer TD-BEM on the `.vol` boundary P1 space
    - first gate: `tests/testVolTdBemConvolutionQuadrature.m` verifies BDF1
      and BDF2 CQ points, positive-real Laplace parameters, residuals,
      real time response, and boundary-data shape checking

10d. Volume-FEM / BEM coupled convolution quadrature, landed 2026-07:
    - `volFemBemCoupledConvolutionQuadrature` is the production-form
      interior/exterior coupled CQ rung. It uses `.vol` tetrahedra as H1/P1
      volume wave FEM and `.vol` boundary triangles as P1 exterior BEM, then
      solves the Johnson-Nedelec / Calderon row at each CQ Laplace point
      `[A+(s/c1)^2 M, -T'Mb; (1/2 Mb-K(s))*T, V(s)] [u_hat; q_hat] = [F_hat; 0]`
    - a volume source pulse drives the interior FEM; the BEM flux density
      radiates to exterior observation points through the matching
      representation `-S(s)q + D(s)Tu`. The inverse CQ FFT recovers interior
      pressure, boundary density, and exterior pressure time histories
    - `K(s)` is the retarded outward-normal principal-value double-layer
      operator, assembled with the same semi-analytic Laplace split as the
      verified frequency-domain Johnson-Nedelec solver. The previous
      `[Mb*T, -V(s)]` single-layer row is still available as
      `CouplingForm="SingleLayerTeaching"` for regression only
    - first gate: `tests/testVolFemBemCoupledConvolutionQuadrature.m` verifies
      an interior-node `.vol`, BDF1/BDF2 positive-real CQ points, coupled
      residuals, real time responses, nonzero interior/exterior fields, and
      nonzero `K(s)` participation

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
