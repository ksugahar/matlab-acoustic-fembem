function sol = fsiCoupledSolve(model, options)
%FSICOUPLEDSOLVE Acoustic fluid-structure interaction (FSI) coupled solve.
%
%   sol = fsiCoupledSolve(model, "Wavenumber", 2.0, ...
%       "LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15);
%   p_s = sol.scatteredAt(points);     % exterior scattered pressure
%   p   = sol.totalAt(points);         % incident + scattered
%
%   % fast exterior for a SPHERE truncation: the exact spherical Helmholtz DtN
%   % / radiating impedance operator instead of the dense Galerkin BEM:
%   sol = fsiCoupledSolve(model, "Wavenumber", 2.0, ..., "ExteriorMethod", "dtn");
%
% The genuine FEM/BEM coupling for acoustics: a solid ELASTIC scatterer
% (vector P1 elasticity FEM in the interior) radiating into an unbounded
% fluid (acoustic BEM in the exterior), a plane wave exp(1i*k*z) incident.
% Unlike a rigid/soft scatterer (which needs no interior PDE - pure BEM),
% an elastic body has internal compressional + shear resonances, so the
% interior displacement field must be solved and coupled to the exterior.
%
% Interface conditions at the boundary Gamma (inviscid fluid):
%   dynamic    : sigma(u).n = -p n      (fluid pressure loads the solid as
%                                        a normal traction, no shear)
%   kinematic  : (1/(rho_f w^2)) dp/dn = u.n   (normal velocities match)
%
% Unknowns (u interior displacement, p_s scattered boundary pressure, q_s
% scattered boundary flux); coupled block system (mirrors the scalar
% Johnson-Nedelec femBemCoupledSolve with a VECTOR interior):
%
%   [ K - w^2 M      G'            0   ] [u  ]   [-G' p_incGamma]
%   [ -rho_f w^2 G   0             Mb  ] [p_s] = [-Mb q_inc     ]
%   [ 0              1/2 Mb - K_k  V_k ] [q_s]   [ 0            ]
%
% with K,M the elasticity stiffness/mass (elasticityMatrices), G the
% coupling int_Gamma mu (n.v), Mb the boundary P1 mass, V_k / K_k the
% acoustic single/double layer. Exterior representation
% p_s = D_k[p_s|Gamma] - S_k[q_s].
%
% EXTERIOR METHOD (options.ExteriorMethod):
%   "bem" (default) - dense Galerkin single/double layer, ANY radiator shape.
%   "dtn"           - the exact spherical Helmholtz DtN (sphericalDtnOperator,
%                     a radiating impedance operator for acoustics, not a
%                     Kelvin boundary). The scattered field is a spherical-harmonic
%                     expansion p_s = Phi c, q_s = Phi diag(Lambda) c, so the
%                     exterior reduces to its (N+1)^2 coefficients c and the
%                     block-3 row becomes the FULL-RANK reduced system
%                       Kdyn u + (G' Phi) c        = -G' p_inc
%                       -rho_f w^2 (Phi' G) u + (Gram diag(Lambda)) c = -Phi' Minc
%                     - no dense N^2 Galerkin assembly (measured ~240x faster
%                     solve on the unit ball). FAIL-LOUD if Gamma is not a sphere.
%
% Speeds/density are RELATIVE TO THE FLUID (fluid c = 1); Lame constants
% mu = DensityRatio*cT^2, lambda = DensityRatio*(cL^2 - 2 cT^2).
%
% VALIDATED (tests/testFsiCoupledSolve): the stiff limit reproduces the rigid
% sphere to ~1e-3 (the formulation gate, independent of interior resolution),
% and the elastic field CONVERGES to elasticSphereScattering under mesh
% refinement (25% coarse -> 7% fine at kR = 2; P1 interior elasticity is the
% accuracy-limiting factor, the coupling is exact). The "dtn" exterior gives
% the SAME field (stiff 3.8e-3 vs rigid, DtN-vs-BEM agreement 7e-4) with the
% dense assembly skipped; the DtN operator itself is exact per multipole
% (2.6e-5 on an independent point-source Dirichlet->Neumann check, degree 10).

arguments
    model (1,1) FemBemModel
    options.Wavenumber (1,1) double {mustBePositive}
    options.LongitudinalSpeed (1,1) double {mustBePositive} = 2.0
    options.ShearSpeed (1,1) double {mustBeNonnegative} = 1.0
    options.DensityRatio (1,1) double {mustBePositive} = 1.5
    options.FluidDensity (1,1) double {mustBePositive} = 1.0
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 7
    options.ExteriorMethod (1,1) string ...
        {mustBeMember(options.ExteriorMethod, ["bem", "dtn"])} = "bem"
    options.DtnDegree (1,1) double {mustBeInteger} = -1
    options.Incident = []          % [] = plane wave; else struct(value, grad)
end

k = options.Wavenumber;
omega = k;                              % fluid c = 1
rhoF = options.FluidDensity;
cL = options.LongitudinalSpeed;
cT = options.ShearSpeed;
rhoS = options.DensityRatio;
mu = rhoS * cT^2;
lamElastic = rhoS * (cL^2 - 2 * cT^2);

mesh = model.mesh;
surface = model.surface;
nV = size(mesh.vtx, 1);
ids = surface.volNodeIds;               % boundary(compact) -> volume node
nB = numel(ids);

% ---- interior: vector elasticity dynamic stiffness ----
[Ks, Ms] = elasticityMatrices(mesh, lamElastic, mu, rhoS);
Kdyn = Ks - omega^2 * Ms;

% ---- incident field (plane wave exp(ikz) by default; a custom Incident struct
% with .value(X) and .grad(X) enables e.g. a phased-array superposition) ----
if isempty(options.Incident)
    incValue = @(X) exp(1i * k * X(:, 3));
    incGrad  = @(X) [zeros(size(X, 1), 2), 1i * k * exp(1i * k * X(:, 3))];
else
    incValue = options.Incident.value;
    incGrad  = options.Incident.grad;
end

% ---- interface coupling + incident data (common to both exterior methods) ----
[Mb, ~] = SurfaceP1Space(surface).mass();
[G, Minc] = interfaceCoupling(surface, ids, nV, incGrad);
pincB = incValue(surface.vtx);

% ---- exterior close + coupled block solve (method-dependent) ----
% co bundles the shared coupling data (interior dynamic stiffness, interface
% coupling G, boundary mass, incident loads) so each exterior method reads as a
% self-contained block assembly + solve.
co = struct("Kdyn", Kdyn, "G", G, "Mb", Mb, "pincB", pincB, "Minc", Minc, ...
    "surface", surface, "nV", nV, "nB", nB, "rhoF", rhoF, "omega", omega);
switch options.ExteriorMethod
    case "bem"
        [u, psG, qs, residual, rhsNorm, dtnInfo] = fsiBemExterior(co, k, options.QuadratureOrder);
    case "dtn"
        [u, psG, qs, residual, rhsNorm, dtnInfo] = fsiDtnExterior(co, k, options.DtnDegree);
end

q = options.QuadratureOrder;
sol = struct();
sol.kind = "acoustic_fluid_structure_interaction_coupled_solve";
sol.policy = "vector_elasticity_fem_interior_" + options.ExteriorMethod + "_exterior";
sol.wavenumber = k;
sol.exteriorMethod = options.ExteriorMethod;
sol.dtn = dtnInfo;
sol.longitudinalSpeed = cL;
sol.shearSpeed = cT;
sol.densityRatio = rhoS;
sol.interiorDisplacement = u;
sol.surfacePressure = psG;
sol.surfaceFlux = qs;
sol.residualNorm = residual;
sol.scatteredAt = @(points) ...
    doubleLayerPotentialMatrix(surface, points, k, q) * psG ...
    - singleLayerPotentialMatrix(surface, points, k, q) * qs;
sol.totalAt = @(points) incValue(points) + sol.scatteredAt(points);
sol.checks = struct( ...
    "solveResidualSmall", residual <= 1e-8 * max(1, rhsNorm), ...
    "fieldComplex", ~isreal(psG), ...
    "exteriorWellPosed", ~dtnInfo.used || dtnInfo.gramCondition < 1e12);
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end


function [G, Minc] = interfaceCoupling(surface, ids, nV, incGrad)
%INTERFACECOUPLING G_ij = int_Gamma mu_i (n . phi_struct_j) and the incident
% normal-flux load Minc_i = int_Gamma mu_i (grad p_inc . n), the incident
% gradient (incGrad: X -> 1x3) dotted with the outward normal at the centroid.
signs = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri;
vtx = surface.vtx;
nB = size(vtx, 1);
massTri = (ones(3) + eye(3)) / 12;      % P1 triangle mass / area

nEntries = 27 * size(tri, 1);
Grow = zeros(nEntries, 1); Gcol = zeros(nEntries, 1); Gval = zeros(nEntries, 1);
Minc = zeros(nB, 1);
cursor = 1;
for t = 1:size(tri, 1)
    lc = tri(t, :);
    X = vtx(lc, :);
    cr = cross(X(2, :) - X(1, :), X(3, :) - X(1, :));
    area = 0.5 * norm(cr);
    nrm = signs(t) * cr / norm(cr);     % outward normal
    vid = ids(lc);                       % volume node ids
    for a = 1:3
        for b = 1:3
            for c = 1:3
                Grow(cursor) = lc(a);
                Gcol(cursor) = 3*vid(b) - 3 + c;
                Gval(cursor) = area * massTri(a, b) * nrm(c);
                cursor = cursor + 1;
            end
        end
    end
    Xc = mean(X, 1);
    gradC = incGrad(Xc);                 % 1x3 incident gradient at the centroid
    qinc = gradC * nrm.';                % (grad p_inc) . n_outward
    Minc(lc) = Minc(lc) + (area / 3) * qinc;
end
G = sparse(Grow, Gcol, Gval, nB, 3 * nV);
end


function [u, psG, qs, residual, rhsNorm, dtnInfo] = fsiBemExterior(co, k, quadOrder)
%FSIBEMEXTERIOR Dense Galerkin single/double layer exterior; nodal (u, p_s, q_s)
% with the exterior Calderon row (1/2 Mb - K_k) p_s + V_k q_s = 0 (any shape).
Vk = GalerkinSingleLayer(co.surface, "Wavenumber", k, "QuadratureOrder", quadOrder).matrix;
Kk = GalerkinDoubleLayer(co.surface, "Wavenumber", k, "QuadratureOrder", quadOrder).matrix;
ZuB = sparse(3 * co.nV, co.nB);
ZBB = sparse(co.nB, co.nB);
lhs = [ co.Kdyn,                       co.G.',           ZuB;
        -co.rhoF*co.omega^2*co.G,      ZBB,              co.Mb;
        ZuB.',                         0.5*co.Mb - Kk,   Vk ];
rhs = [ -co.G.' * co.pincB; -co.Minc; zeros(co.nB, 1) ];
x = lhs \ rhs;
u   = x(1:3 * co.nV);
psG = x(3 * co.nV + (1:co.nB));
qs  = x(3 * co.nV + co.nB + (1:co.nB));
residual = norm(lhs * x - rhs);
rhsNorm = norm(rhs);
dtnInfo = struct("used", false);
end


function [u, psG, qs, residual, rhsNorm, dtnInfo] = fsiDtnExterior(co, k, dtnDegree)
%FSIDTNEXTERIOR Exact spherical Helmholtz DtN exterior (radiating impedance).
% The scattered field is a spherical-harmonic expansion p_s = Phi c,
% q_s = dp_s/dn = Phi diag(Lambda) c; reducing the exterior to its nModes
% coefficients c keeps the system FULL RANK and replaces the dense N^2 Galerkin
% assembly by a low-rank surface operator.  Fail-loud if Gamma is not a sphere.
%   Kdyn u               + (G' Phi) c              = -G' p_inc
%   -rhoF w^2 (Phi' G) u + (Gram diag(Lambda)) c   = -Phi' Minc
dtn = sphericalDtnOperator(co.surface, "Wavenumber", k, "Degree", dtnDegree);
Phi = dtn.harmonics;                 % nB x nModes
lam = dtn.modeEigenvalues;           % nModes x 1 (Lambda per column)
nM = dtn.numModes;
lhs = [ co.Kdyn,                             co.G.' * Phi;
        -co.rhoF*co.omega^2*(Phi.' * co.G),  dtn.gram .* lam.' ];
rhs = [ -co.G.' * co.pincB; -(Phi.' * co.Minc) ];
x = lhs \ rhs;
u = x(1:3 * co.nV);
c = x(3 * co.nV + (1:nM));
psG = Phi * c;                       % nodal scattered pressure
qs  = Phi * (lam .* c);              % nodal scattered flux T[p_s]
residual = norm(lhs * x - rhs);
rhsNorm = norm(rhs);
dtnInfo = struct("used", true, "degree", dtn.degree, "numModes", dtn.numModes, ...
    "radius", dtn.radius, "sphericityDeviation", dtn.sphericityDeviation, ...
    "gramCondition", dtn.gramCondition);
end
