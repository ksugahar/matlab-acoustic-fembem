function tests = testBemLayerSelfConsistency
%TESTBEMLAYERSELFCONSISTENCY Single- + double-layer self-consistency (Calderon / Green).
%
% Method of manufactured solutions, NO analytic partial-wave series: put a point
% source at x0 INSIDE the unit sphere, so u = G(.,x0) is an exact radiating
% solution in the EXTERIOR.  Its exterior Cauchy data (p = u|_S Dirichlet,
% q = du/dn|_S Neumann, n outward) must satisfy, to P1 discretization error, two
% identities that couple the single (V, S) and double (K, D) layer operators:
%
%   (A) boundary Calderon identity   V q = (K - 1/2 M) p     -- M appears EXPLICITLY
%   (B) Green representation         u(x) = (D p - S q)(x)   exterior; 0 interior
%
% (A) is a mass-matrix consistency guard: dropping the boundary P1 mass M breaks
% it by ~45%.  This is the operator-suite analogue of the CQ single-layer RHS
% fix (which needed V q = M ghat, not the raw nodal ghat).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function [surf, p, q, Gfun, k, s, c0] = manufacturedExteriorField()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_fine.vol");
surf = VolMesh(volFile).boundary();
c0 = 1.0; k = 1.8; s = -1i * c0 * k;             % kernel exp(+i k r)/(4 pi r)

Y = surf.vtx;
nrm = Y ./ vecnorm(Y, 2, 2);                     % outward unit normal (unit sphere)
x0 = [0.2 0 0];                                  % point source INSIDE the sphere
Gfun = @(X) exp(1i*k*vecnorm(X-x0,2,2)) ./ (4*pi*vecnorm(X-x0,2,2));
r = vecnorm(Y - x0, 2, 2);
p = Gfun(Y);                                     % Dirichlet trace
dirdot = sum((Y - x0).*nrm, 2) ./ r;             % ((Y-x0)/r . n)
q = p .* (1i*k - 1./r) .* dirdot;                % Neumann trace (outward n)
end


function testExteriorCalderonIdentityNeedsMass(testCase)
% V q = (K - 1/2 M) p holds to discretization error; the SAME identity with the
% boundary mass M dropped is grossly violated (proves M is structurally required).
[surf, p, q, ~, ~, s, c0] = manufacturedExteriorField();
order = 7;
V = laplaceSingleLayerGalerkin(surf, s, c0, order);
K = laplaceDoubleLayerGalerkin(surf, s, c0, order);
[M, ~] = SurfaceP1Space(surf).mass();

den = norm(V * q);
relWithMass = norm(V*q - (K - 0.5*M)*p) / den;
relNoMass   = min(norm(V*q - K*p), norm(V*q + K*p)) / den;

verifyLessThan(testCase, relWithMass, 1.5e-2);   % measured 5.5e-3 (Calderon holds)
verifyGreaterThan(testCase, relNoMass, 1.0e-1);  % measured 0.45 (drop M -> breaks)
end


function testGreenRepresentationReconstructsExteriorAndVanishesInterior(testCase)
% u(x) = (D p - S q)(x): reproduces the manufactured field at exterior points,
% and vanishes at interior points (the exterior representation's null-field).
[surf, p, q, Gfun, ~, s, c0] = manufacturedExteriorField();
order = 7;
extObs = [0 0 3; 2.5 0 0; 0 3 0];
intObs = [0 0 0.2; -0.3 0 0; 0 0.25 0];

Sext = laplaceSingleLayerPotential(surf, extObs, s, c0, order);
Dext = laplaceDoubleLayerPotential(surf, extObs, s, c0, order);
Sint = laplaceSingleLayerPotential(surf, intObs, s, c0, order);
Dint = laplaceDoubleLayerPotential(surf, intObs, s, c0, order);
uext = Gfun(extObs);
uint = Gfun(intObs);

relExterior = norm(Dext*p - Sext*q - uext) / norm(uext);
relInterior = norm(Dint*p - Sint*q) / norm(uint);
verifyLessThan(testCase, relExterior, 5e-2);     % measured 2.2e-2 (representation)
verifyLessThan(testCase, relInterior, 3e-2);     % measured 3.4e-3 (null-field)
end
