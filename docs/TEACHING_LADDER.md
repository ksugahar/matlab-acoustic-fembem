# FEM/BEM Teaching Ladder

This document preserves the detailed solver narrative for
`matlab-acoustic-fembem`: the readable MATLAB API, the cross-validation ladder,
and the acoustic/FSI examples that are too long for the top-level README.

This repository is an open educational MATLAB codebase. It is not a production
solver and it is not meant to compete with NGSolve on speed. Its job is to make
the mathematics readable and reproducible:

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

## License and Gypsilab Relationship

This repository is released under the GNU GPL v3.0.  Gypsilab itself is
GPL-3.0 software, and this project intentionally follows its readable MATLAB
style.  The repository does not vendor the Gypsilab source tree; optional
cross-checks against a separately installed Gypsilab checkout are kept behind
explicit verification helpers.  Please cite Gypsilab and its authors when this
teaching code helps your FEM/BEM work.

## Basic API

```matlab
repoRoot = pwd;                 % or the cloned repository path
addpath(genpath(fullfile(repoRoot, "matlab_api")));

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
quadrature. Cross-checked against TWO independent codes on the same sphere
meshes:

- the real Gypsilab (`validation/verifyGalerkinAgainstGypsilab.m`): operator
  1.1e-4, capacitance 1.4e-5 at gss 3 - tight because both codes use the
  same 3-point test quadrature, so that error cancels;
- NGSolve's `ngsolve.bem` (`validation/verifyGalerkinAgainstNgsolve.m`,
  reading committed reference matrices from
  `validation/exportNgsolveBemReference.py`): V 3.8e-4 / K 3.4e-3 at gss 7,
  capacitance 1.6e-4, mass identical to 1e-16. The reference is
  Sauter-Schwab at intorder 16 (self-converged to 1e-8), so this one
  measures the TRUE assembly error - the gss 1/3/7 sweep converges
  7e-2 / 7e-3 / 4e-4 against it. Conventions match exactly (same
  outward-normal principal-value K: K[1] = -1/2 to 1e-9). This check runs
  in the test suite (`testNgsolveBemCrossCheck`) from the committed .mat
  artifacts alone.

The final rung (stage 5) is the Johnson-Nedelec coupled FEM/BEM
open-boundary solve: P1 FEM inside the meshed volume, Galerkin BEM for the
unbounded exterior, continuous trace and flux on the boundary:

```matlab
m = FemBemModel("fixtures/mesh_topology/unit_sphere_fine.vol");
sol = femBemCoupledSolve(m);        % -div(c grad u) = f inside, Laplace outside
sol.u                               % interior nodal solution
sol.lambda                          % exterior normal flux on the boundary
u = sol.exteriorPotentialAt([3 0 0]);
```

The coupled system is the classic pair `A u - T' M lambda = F` and
`(1/2 M - K) T u + V lambda = 0`, with `GalerkinDoubleLayer` supplying the
outward-normal principal-value K (sphere gates: K[1] = -1/2 exact to six
digits, K[Y_1] = -1/6 to 0.3%). The unit-ball source f = 1 reproduces the
analytic radial solution u = 1/2 - r^2/6 (trace mean -2.3% on the fine
fixture, flux conservation to 1e-3, exterior potential 3.5%; all
geometry-faceting dominated).

## Vector Coupling Teaching Path (H(curl)/RWG, ladder stage 6)

The Maxwell-precursor layer makes the vector trace coupling concrete and
locks it to machine precision:

```matlab
m = FemBemModel("mesh.vol");
C = m.rwg.rotatedTraceMap(m.hcurl);   % n x u|_Gamma in RWG coords: c = C*alpha
[G, d] = m.rwg.gram();                % RWG mass, exact with the 3-point rule
L = RwgSingleLayer(m.rwg);            % static EFIE / partial-inductance kernel
A = L.vectorPotentialAt(c, points);   % analytic vector potential of sum c f
```

The FEEC identity behind `rotatedTraceMap` is tested pointwise on every
boundary edge: the rotated tangential trace of the volume Nedelec0 edge
function IS the RWG function of the same edge,

```
n x (gamma_t N_E) = -signOut * triEdgeSigns * sigma_pm * f_e / l_e
```

so a volume H(curl) field crosses to the boundary RWG space exactly, with
no interpolation. `RwgSingleLayer` reuses the verified P1 panel integrals
(RWG components are affine per triangle - no new singular math) and is
validated by the uniformly magnetized sphere: K = z_hat x n reproduces
A = (1/3) z_hat x x inside and the dipole field outside (coarse 3.1%,
fine 1.3%). The full vector transmission solve (eddy-current FEM/BEM with
the vector Calderon operators) is intentionally left to NGSolve.BEM.

## Acoustic Helmholtz Rung (3-way validated)

The exterior single-layer solve extends to acoustics through the same
low-frequency-stable kernel split (analytic Laplace panels + smooth
`(exp(1i*k*r) - 1)/(4*pi*r)` correction; `e^{+ikr}` convention carries the
radiation condition):

```matlab
surface = VolMesh("fixtures/mesh_topology/unit_sphere_fine.vol").boundary();
g = -exp(1i * k * surface.vtx(:, 3));        % sound-soft: p_scat = -p_inc
sol = singleLayerDirichletSolve(surface, g, "Wavenumber", k);
p = sol.potentialAt(points);                 % scattered field outside
```

Validated THREE ways on the same meshes
(`validation/verifyHelmholtzAgainstNgsolve.m`,
`tests/testHelmholtzScattering.m`):

- analytic: an interior point source (EXACT gate - no truncation, pure
  discretization error; `acousticPointSource`) and the sound-soft sphere
  partial-wave series (`softSphereScattering`; truncation tail reported);
- NGSolve's `ngsolve.bem` Helmholtz references: committed .mat with the
  complex V, its intorder self-convergence, and NGSolve's own solved
  GetPotential probe values (`validation/exportNgsolveBemHelmholtzReference.py`).

The reading that matters: the two CODES agree 10-30x tighter with each
other (probe cross-code 3e-4..6e-3) than either agrees with the true
sphere (1-10%, faceting-dominated, improving coarse -> fine), proving the
analytic deviation is geometry, not implementation. The constant-mode
eigenvalue `lambda_0 = sin(k) e^{+ik} / k` pins the time convention on
both sides. The first-kind V_k equation is taught WITH its
irregular-frequency caveat (unit sphere: first interior Dirichlet
eigenvalue at kR = pi); CHIEF / Burton-Miller is the next acoustic rung.

### Sonic-Crystal Chain (multi-body, 4-leg validated)

Multiple scatterers need no new machinery - the all-pairs Galerkin
assembly handles any number of closed surfaces in one `.vol`. The
teaching fixture is a five-sphere sound-soft chain (R = 0.3, lattice
constant d = 1.5; `validation/makeSoftSphereChainFixture.py`):

```matlab
surface = VolMesh("fixtures/mesh_topology/soft_sphere_chain_5.vol").boundary();
sol = singleLayerDirichletSolve(surface, -exp(1i*k*surface.vtx(:,3)), ...
    "Wavenumber", k, "QuadratureOrder", 3);
```

Validated FOUR ways (`validation/verifySonicCrystalChain.m`,
`tests/testSonicCrystalChain.m`): the exact interior-point-source gate
(source inside the MIDDLE sphere, 1.8-2.0e-2), Foldy monopole multiple
scattering (`foldyPointScattering`, the low-k analytic-class reference:
2.7e-2 at k = 1.0, honestly degrading to ~15% at kd = pi as the neglected
l >= 1 terms turn on), the ngsolve.bem reference (probe cross-code
4e-4..1.1e-3), and an NGSolve VOLUME-FEM leg (order 2, Dirichlet spheres,
first-order Sommerfeld ABC - a fully different discretization, agreeing
to 4.4-7.5e-2 at the probes).

PHYSICS FINDING (locked as a negative-result test): the sparse free-space
chain shows broadband sub-wavelength attenuation (~3.5-3.9 dB insertion
loss, flat over k = 0.6..3.0) and NO Bragg stop band at k d = pi - all
four legs agree. A duct-confined sonic-crystal model
requires duct confinement / transverse periodicity; the Bloch unit cell +
duct transmission FEM is the declared next rung.

### Rigid Scattering, Irregular Frequencies, CHIEF

Sound-hard scattering uses the total-field double-layer equation
`(1/2 M - K_k) t = M g_inc` (`rigidScatteringSolve`; the representation
loses its single-layer term because dp/dn = 0, so the exterior field is
just `u_inc + D_k[t]`). Cross-checked at the operator level (our K_k vs
ngsolve.bem: 1.8-3.1e-3, conjugate 0.12-1.4) and at the solve level -
the second-kind cross-code agreement (trace 4e-6..1e-4) is 1-2 orders
TIGHTER than the first-kind V solves, the conditioning lesson made
visible. Both codes sit 3.8-13% from the rigid-sphere series (shared
faceting).

The taught trap, locked as a test: at kR = pi (first interior Dirichlet
eigenvalue) the equation is singular - the DISCRETE condition number
stays benign (~29, the faceted eigenvalue shifts) while the answer is
~100% wrong; only the analytic gate sees it. CHIEF (Schenck 1968:
interior null-field rows at jittered interior points, least squares)
restores regular-class accuracy (0.96 -> 0.080); Burton-Miller is the
production alternative (hypersingular operator - intentionally beyond
this lane, see NGSolve.BEM).

### Duct Band Gap (stage 9: the confinement experiment)

The stage-8 finding said free space cannot make the gap; the duct does.
`validation/exportNgsolveDuctBandReference.py` (NGSolve FEM, committed
.mat, locked by `tests/testDuctBandGap.m`) puts the SAME sound-soft
sphere family (R = 0.3, d = 1.5) in a rigid duct (a = 1.0, single-mode
below pi):

- Bloch bands of the unit cell (Floquet phase, order 2; empty-lattice
  analytic gate 2.0e-4): band 1 [2.31, 2.52], gap [2.52, 3.65];
- transmission through the finite 5-cell crystal (one-way ports; empty
  duct transparent to 1e-4): stop (T ~ 1e-4 below band 1, soft
  inclusions cut off long wavelengths) / pass (T up to 0.27 inside
  band 1) / stop (T ~ 1e-3 in the Bragg gap) - aligned with the bands,
  contrast ratio > 100.

Free space (stage 8, 4-leg): no gap. Duct (stage 9): full band
structure. That contrast is the sonic-crystal teaching story.

## Coupled Acoustic FEM/BEM (pure MATLAB, 3 gated cases)

The Johnson-Nedelec solve extends to the acoustic TRANSMISSION problem:
interior P1 Helmholtz FEM (medium 1: `k1`, `rho1`) coupled to the exterior
P1 Galerkin BEM (medium 0) - with NO absorbing boundary anywhere, because
the BEM row IS the exact radiation condition (contrast the stage-8 volume
FEM leg, whose first-order ABC costs 4-7% at the probes):

```matlab
m = FemBemModel("fixtures/mesh_topology/unit_ball_maxh018.vol");
sol = femBemCoupledSolve(m, "Wavenumber", 2.0, ...
    "InteriorWavenumber", 2/0.7, "DensityRatio", 1.2, ...
    "VolumeSource", 0, "IncidentAmplitude", 1);
sol.u                               % interior total pressure (complex)
u_s = sol.exteriorPotentialAt(x);   % scattered field, -S_k + D_k
```

The one missing operator was the Helmholtz double layer:
`K_k = K_0` (the verified analytic panels) `+` smooth correction
`base*(exp(z)(1-z)-1)` through `HelmholtzKernel` SourceNormals - no new
singular math. Sphere spectral gates:
`K_k[Y_l] = 1/2 + 1i k^2 j_l(k) h_l'(k)` to 3-5e-3 at k = 0.5 and
2-3e-2 at k = 2 (faceting class); the k -> 0 limit is the Laplace K to
5e-28.

Three gated cases (`tests/testFemBemHelmholtzCoupling.m`):

1. **k -> 0 regression** - matches the verified Laplace coupled solve to
   1e-9 (density 1e-12).
2. **Acoustic invisibility** (k1 = k0, rho1 = rho0): the exact null gate -
   interior == incident plane wave, scattered == 0, up to discretization:
   4.1e-2 / 2.6e-2 (unit_sphere_fine) -> 1.3e-2 / 7.8e-3 (unit_ball_maxh018),
   locked as a CONVERGENCE assertion, not a loose band.
3. **Anderson (1950) fluid sphere** (c1/c0 = 0.7, rho1/rho0 = 1.2) against
   the partial-wave transmission series (`fluidSphereScattering`,
   log-derivative-stable mode solve): interior 13% -> 4.4%, exterior
   probes 22% -> 7.3% under refinement - the P1 (k1 h)^2 resolution class
   measured, not hidden.

### Frequency-Sweep Time Response

`volFemBemIfftResponse` is the first executable time-domain bridge built on
the real `.vol` FEM/BEM stack:

```matlab
td = volFemBemIfftResponse("fixtures/mesh_topology/unit_tetra.vol", ...
    "NumTime", 16, "FinalTime", 24);
```

It reads the Netgen `.vol` mesh, assembles H1/P1 tetrahedral volume FEM and
boundary P1 Galerkin BEM through `FemBemModel`, solves `femBemCoupledSolve`
on many positive frequency bins, multiplies the transfer function by a real
pulse spectrum, mirrors the spectrum with Hermitian symmetry, and obtains a
real pressure time history by `ifft`.

This is not a sine-wave animation and not yet convolution-quadrature
time-domain BEM. It is the readable frequency-domain FEM/BEM plus inverse-FFT
route toward the full P1 volume-FEM/P1 surface-BEM transient solver.

### Convolution-Quadrature TD-BEM

`volTdBemConvolutionQuadrature` is the first Lubich CQ time-domain BEM rung:

```matlab
td = volTdBemConvolutionQuadrature("fixtures/mesh_topology/unit_tetra.vol", ...
    "NumTime", 8, "TimeStep", 0.6, "Method", "BDF1");
```

It uses the `.vol` boundary triangles as a genuine 3D P1 Galerkin BEM surface,
samples the BDF generating function `delta(zeta)`, evaluates Laplace-domain
retarded single-layer matrices `V(s)` at `s = delta(zeta)/dt` with `Re(s)>0`,
solves `V(s) qhat = ghat`, and FFT-recovers both the boundary density and
exterior pressure time histories.

This is different from the drum GIF rung. The current drum visualizer is a
3D axisymmetric physical picture rendered on an r-z slice: disk and side
sources are integrated over azimuth, but the mesh is not yet a 3D `.vol`
P1 volume-FEM/P1 surface-BEM drum. A true drum next step is to generate a
cylindrical drum `.vol`, then feed its boundary triangles to this CQ-BEM lane
and later add the coupled volume-FEM interior equation.

`volFemBemCoupledConvolutionQuadrature` adds that first volume-FEM interior:

```matlab
td = volFemBemCoupledConvolutionQuadrature( ...
    "fixtures/mesh_topology/four_tet_interior_node.vol", ...
    "NumTime", 8, "TimeStep", 0.6);
```

At each CQ Laplace point it solves the readable block system

```text
[ A+(s/c1)^2 M       -T' Mb ] [ u_hat ] = [ F_hat ]
[ (1/2 Mb-K(s)) T     V(s)  ] [ q_hat ]   [   0   ]
```

The volume source pulse drives the interior P1 FEM, the boundary trace is
coupled to the exterior P1 BEM flux density, and the exterior pressure is
recovered at observation points by `-S(s)q + D(s)Tu`. This is the
Calderon/Johnson-Nedelec coupled CQ system with the retarded double-layer
`K(s)`. The earlier `[Mb*T, -V(s)]` single-layer row remains available only as
`CouplingForm="SingleLayerTeaching"` for regression and classroom contrast.

## Adjoint Automatic Differentiation (wavefront synthesis)

Reverse-mode AD *through* the BEM solve, the readable way: differentiate
the equation, not the LU. `acousticFocusAdjoint` designs a phased array's
complex amplitudes (phases) to focus acoustic energy at a target point
behind a rigid scatterer - the gradient of the focused intensity
`J = |u(target)|^2` comes from ONE adjoint (transpose) solve, independent
of the number of sources:

```matlab
surface = VolMesh("fixtures/mesh_topology/unit_sphere_fine.vol").boundary();
sources = ...;                 % phased array, one complex amplitude each
res = acousticFocusAdjoint(surface, sources, [0 0 2.5], 2.0, amplitudes, ...
    "GradientCheck", true);
res.gradientReal / res.gradientImag   % dJ/dRe(p), dJ/dIm(p) - one solve
res.ascentDirection                    % steepest ascent = 2 u conj(w)
res.gradientCheckRelError              % vs central finite differences
```

Because the field is affine in the amplitudes, `u = w * p` with the
sensitivity row `w = S0 + lambda' M S` assembled from the adjoint solve
`A' lambda = d0'`. Measured (`tests/testAcousticFocusAdjoint.m`): forward
affine residual 4e-18 (w is exact), adjoint gradient vs central finite
differences 1.7e-10, and a gradient ascent monotonically focuses the
energy (the wavefront-synthesis proof). The objective is non-holomorphic
(`|.|^2`), so the ascent direction is the Wirtinger derivative
`dJ/dconj(p) = 2 u conj(w)` - NOT `2 conj(u) w` (a sign trap that a
near-zero starting field hides and a nonzero start exposes).

This is the seed of acoustic radiation-force / thrust design: the
intensity objective becomes the net radiation force (below), and the same
adjoint gives dF/d(phases) for wavefront-synthesised propulsion.

## Acoustic Radiation Force (ultrasonic thrust / manipulation)

`acousticRadiationForce` computes the net time-averaged force on a
scatterer - the SECOND-order functional of the linear field, as the
closed-surface integral of the Brillouin radiation-stress tensor over a
control sphere:

```matlab
sol = rigidScatteringSolve(surface, "Wavenumber", 2.0, "QuadratureOrder", 7);
pField = @(X) exp(1i*2*X(:,3)) + doubleLayerPotentialMatrix(surface, X, 2, 7) * sol.trace;
res = acousticRadiationForce(pField, 2.0);
res.force            % 1x3 net force; res.forceFunction = Y_p
```

The same post-processor takes the analytic partial-wave series OR the BEM
total field (geometry-general). Gated the formula-free way
(`tests/testAcousticRadiationForce.m`): the force is INDEPENDENT of the
control radius to ~1e-10 (the tensor is divergence-free in the source-free
fluid - the primary correctness gate, no external formula needed), pushes
downstream (F_z > 0), is axisymmetric, and gives `Y_p(kR=2) = 0.752`
(King 1934); the BEM force matches the series to ~5% (faceting).

Physical, air @ 40 kHz (c = 343 m/s, rho = 1.2, lambda = 8.58 mm): the
dimensionless kR = 2 is a 2.73 mm sphere, and
`F_z = Y_p * pi * Rs^2 * <E>` with `<E> = p_rms^2/(rho c^2)`:

| SPL | p (Pa) | force on Rs = 2.73 mm |
|-----|--------|-----------------------|
| 120 dB | 28 | 0.05 uN |
| 140 dB | 283 | 5 uN |
| 160 dB | 2828 | 0.5 mN |

That is the acoustic-levitation / micro-manipulation regime (a mm bead's
weight is ~10 uN, reached near 150 dB - consistent with real 40 kHz
levitator arrays). Feed this force as the objective to the adjoint above
for dF/d(phases) - wavefront-synthesised thrust design.

## Elastic Sphere / Fluid-Structure Reference (FSI)

`elasticSphereScattering` is the Faran (1951) solid ELASTIC sphere - the
analytic reference for the acoustic fluid-structure INTERACTION case,
where FEM/BEM coupling truly earns its keep (an elastic scatterer has
internal compressional + shear resonances that rigid/soft spheres miss):

```matlab
ref = elasticSphereScattering(k, 1.0, points, ...
    "LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15);
Y = acousticRadiationForce(@(X) elasticSphereScattering(k,1,X,...).total, k);
```

Validated (`tests/testElasticSphereScattering.m`) against TWO independent
references: ShearSpeed -> 0 reproduces the Anderson fluid sphere to
1.3e-10, a stiff/heavy solid reproduces the rigid sphere. A sign bug in
the radial-stress boundary term passed the stiff->rigid limit but failed
the fluid limit by 20-70% - only the independent Anderson reference
exposed it (gate on both). The elastic radiation force shows the
resonance: a lucite-like sphere has Y_p(kR) peaking at 2.83 near kR = 3
(well above the rigid ~0.75), the internal resonance rigid/soft cannot
produce. ## Acoustic FSI Coupled Solve (fluid-structure interaction)

`fsiCoupledSolve` is the genuine acoustic FEM/BEM coupling: a solid ELASTIC
scatterer (vector P1 elasticity FEM interior) radiating into an unbounded
fluid (acoustic BEM exterior), coupled through the interface conditions
(fluid pressure loads the solid as a normal traction; fluid and solid
normal velocities match). This is the case where FEM/BEM coupling is FOR -
the interior displacement field (internal compressional + shear resonances)
must be solved, unlike a rigid/soft scatterer (pure BEM).

```matlab
m = FemBemModel("fixtures/mesh_topology/unit_ball_maxh018.vol");
sol = fsiCoupledSolve(m, "Wavenumber", 2.0, ...
    "LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15);
p = sol.totalAt(points);           % exterior total pressure
u = sol.interiorDisplacement;      % solved elastic field inside
```

The coupled block system mirrors the scalar Johnson-Nedelec
femBemCoupledSolve with a VECTOR interior (unknowns interior displacement u,
scattered boundary pressure p_s, scattered flux q_s):

```
[ K - w^2 M      G'            0   ] [u  ]   [-G' p_inc]
[ -rho_f w^2 G   0             Mb  ] [p_s] = [-Mb q_inc]
[ 0              1/2 Mb - K_k  V_k ] [q_s]   [ 0       ]
```

Validated (`tests/testFsiCoupledSolve.m`) against the analytic elastic
sphere: the stiff limit reproduces the rigid sphere to ~1e-3 (the
formulation gate - the interface + BEM coupling is exact, independent of
interior resolution), and the elastic field CONVERGES to
elasticSphereScattering under mesh refinement (25% coarse -> 7% fine at
kR = 2; the P1 interior elasticity - especially the shear field at
kR ~ wavelength - is the accuracy-limiting factor, not the coupling).

This is the acoustic analogue of the lab's core EM coupling problems
(nonlinear iron + open boundary, eddy currents in a conductor + open
boundary), where FEM handles the complex bounded interior and BEM (or the
BEM/Sommerfeld/DtN/high-order-impedance lane) handles the acoustic open
boundary.  In the acoustic lane, do not call the boundary Kelvin.

### Spherical DtN Fast Path (Radiating-Impedance Operator on the Sphere)

When the acoustic truncation surface is a SPHERE, the exterior needs no dense
Galerkin BEM at all: the exact spherical Helmholtz Dirichlet-to-Neumann map
(`sphericalDtnOperator`) imposes the radiation condition as a low-rank
surface radiating-impedance operator.

```matlab
sol = fsiCoupledSolve(model, "Wavenumber", 2.0, ..., "ExteriorMethod", "dtn");
```

The scattered field is a spherical-harmonic expansion `p_s = Phi c`
(`Lambda_n = k h_n'(kR)/h_n(kR)`), so the exterior reduces to its `(N+1)^2`
coefficients, the coupled system stays FULL RANK (a nodal low-rank DtN would
leave `p_s` under-determined in the mesh null space), and the singular `N^2`
Galerkin assembly is skipped entirely. Measured on the unit ball (kR = 2):

- the DtN operator is EXACT per multipole - an INDEPENDENT point-source
  Dirichlet->Neumann check gives 2.6e-5 at degree 10, converging
  7e-3 / 1e-3 / 2e-4 / 3e-5 at degree 4 / 6 / 8 / 10 to the ~1e-5 mesh floor;
- the coupled FSI field matches the rigid sphere in the stiff limit to
  3.8e-3 (dense-BEM leg 4.2e-3; the DtN and BEM legs agree to 7e-4), and
  converges to the elastic sphere in the same P1-interior band;
- the dense order-7 Galerkin assembly (~85 s on the unit ball) is gone - the
  DtN solve is sub-second (two orders of magnitude faster).

This is the acoustic "sphere where spherical DtN applies" rule: a spherical
truncation gets the fast exact operator; a general radiator stays on the
dense BEM or uses a high-order impedance/absorbing boundary. The path is
FAIL-LOUD - `sphericalDtnOperator` raises on a non-spherical surface rather
than silently falling back to the sphere.

### Elastic-Bead Thrust Design (wavefront synthesis through the FSI)

`elasticThrustAdjoint` closes the ultrasonic-thrust story: design a phased
array's complex amplitudes to steer the acoustic RADIATION FORCE on a solid
ELASTIC bead - through the full FSI solve, because an elastic bead has internal
compressional + shear resonances a rigid/soft bead cannot show.

```matlab
res = elasticThrustAdjoint(model, sources, k, amplitudes, "GradientCheck", true);
res.force            % 1x3 net radiation force on the elastic bead
res.ascentDirection  % dF_z/dconj(p): raise the thrust
res.forceForm        % {Q1,Q2,Q3}: reuse in a design loop, F_i(p) = p^H Q_i p
```

Three readable steps compose the earlier rungs (adjoint AD + radiation force +
FSI):

1. **The field is linear in the amplitudes** - one FSI solve per source
   (incident = a unit monopole), then `p_total(X) = fieldAt(X) * p` with the
   double/single-layer potentials SHARED over sources.
2. **The force is a quadratic form** - the radiation force is a second-order
   functional of the linear field, so each component is `F_i(p) = p^H Q_i p`
   (Q_i Hermitian), assembled from the Brillouin stress on the control sphere.
3. **The gradient is one matvec** - the Wirtinger steepest ascent of the real
   force is `dF_i/dconj(p) = 2 Q_i p`. Force AND its full gradient are closed
   forms once Q_i is built; a design loop REUSES `res.forceForm` with NO
   re-solve.

Validated (`tests/testElasticThrustAdjoint.m`): the quadratic form reproduces a
vectorised direct Brillouin integral to ~1e-15 (an independent assembly -
outer-product form vs scalar-field integral, which caught a `conj(z) z.'` vs
`z z'` conjugation error), that matches the golden `acousticRadiationForce`, the
force is control-radius independent (div T = 0), and the Wirtinger gradient
matches central finite differences to ~1e-13; a step along the ascent raises
F_z. The N FSI solves are the only real cost (sub-second each via the DtN
exterior), so build the reusable `forceForm` once (~15 s) and optimise the
phases freely.

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
run(fullfile(repoRoot, "run_tests.m"))
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
- the validation log is recorded under a local `_crossval` directory

The test suite checks that the catalog has exactly 100 unique cases, that all
verified cases have example scripts, and that validation logs are recorded.

Current progress:

- `100 / 100` example scripts present
- `100 / 100` verified
- verified groups: all 10 categories
- validation log:
  `_crossval/gypsilab_mesh_topology_10of100_20260624.md`
  `_crossval/gypsilab_remaining_90of100_20260624.md`
  `_crossval/gypsilab_class_api_refactor_20260703.md`

For acoustic FEM/BEM cases, the public gates rely on analytic references,
NGSolve/NGSolve.BEM reference artifacts, and committed MATLAB fixtures.
