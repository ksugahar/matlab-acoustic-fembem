function tests = testFemBemCoupledSolve
%TESTFEMBEMCOUPLEDSOLVE Johnson-Nedelec coupled solve, ladder stage 5.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDoubleLayerSphereSpectralIdentities(testCase)
% unit sphere, outward-normal kernel: K[1] = -1/2, K[Y_1] = -1/6.
surface = fixtureSurface("unit_sphere_coarse.vol");
K = GalerkinDoubleLayer(surface);
[M, ~] = SurfaceP1Space(surface).mass();
one = ones(size(surface.vtx, 1), 1);
z = surface.vtx(:, 3);

k0 = (one.' * (K.matrix * one)) / (one.' * (M * one));
k1 = (z.' * (K.matrix * z)) / (z.' * (M * z));
verifyEqual(testCase, k0, -0.5, "AbsTol", 5e-3);
verifyEqual(testCase, k1, -1/6, "AbsTol", 5e-3);
end


function testExteriorBieModesAreConsistent(testCase)
% direct exterior BIE V*lambda = (K - 1/2 M)*u for u_e = 1/r and z/r^3.
surface = fixtureSurface("unit_sphere_coarse.vol");
V = GalerkinSingleLayer(surface);
K = GalerkinDoubleLayer(surface);
[M, ~] = SurfaceP1Space(surface).mass();
one = ones(size(surface.vtx, 1), 1);
z = surface.vtx(:, 3);

res0 = norm(V.matrix * (-one) - (K.matrix - 0.5 * M) * one);
verifyLessThan(testCase, res0, 0.04 * norm((K.matrix - 0.5 * M) * one));
res1 = norm(V.matrix * (-2 * z) - (K.matrix - 0.5 * M) * z);
verifyLessThan(testCase, res1, 0.05 * norm((K.matrix - 0.5 * M) * z));
end


function testUnitBallSourceMatchesAnalyticSolution(testCase)
% -Delta u = 1 in the unit ball, Laplace outside, u -> 0 at infinity:
%   u(r) = 1/2 - r^2/6,  u_Gamma = 1/3,  lambda = -1/3.
% Bands locked from the 2026-07-03 measurement on unit_sphere_fine
% (18 interior nodes; errors are geometry-faceting dominated).
model = FemBemModel(fixturePath("unit_sphere_fine.vol"));
verifyGreaterThan(testCase, ...
    size(model.mesh.vtx, 1) - numel(model.mesh.traceNodeIds), 0, ...
    "fixture must have interior FEM nodes");

sol = femBemCoupledSolve(model);

verifyEqual(testCase, sol.status, "ok");
verifyEqual(testCase, sol.kind, "johnson_nedelec_coupled_fem_bem_solve");
verifyLessThan(testCase, sol.residualNorm, 1e-10);

r = sqrt(sum(model.mesh.vtx.^2, 2));
uExact = 0.5 - r.^2 / 6;
verifyLessThan(testCase, max(abs(sol.u - uExact)), 0.03);

[M, ~] = SurfaceP1Space(model.surface).mass();
area = sum(model.surface.areas());
traceMean = sum(M * sol.trace) / area;
lambdaMean = sol.totalExteriorFlux / area;
verifyEqual(testCase, traceMean, 1/3, "AbsTol", 0.02);
verifyEqual(testCase, lambdaMean, -1/3, "AbsTol", 0.02);

% discrete conservation: exterior flux balances the volume source
verifyEqual(testCase, sol.totalExteriorFlux, -sol.volumeSourceTotal, ...
    "RelTol", 1e-3);
end


function testExteriorPotentialMatchesFarField(testCase)
% outside the ball the exact solution is (1/3)/r.
model = FemBemModel(fixturePath("unit_sphere_fine.vol"));
sol = femBemCoupledSolve(model);
points = [3 0 0; 0 0 5];
u = sol.exteriorPotentialAt(points);
uExact = (1/3) ./ [3; 5];
verifyEqual(testCase, u, uExact, "RelTol", 0.06);
end


function testCoefficientScalesTheInteriorSolution(testCase)
% with c = 2 the interior solution is u = 1/3 + (1/2 - r^2/6 - 1/3)/2
% (exterior part unchanged: lambda and u_Gamma stay -1/3 and 1/3).
model = FemBemModel(fixturePath("unit_sphere_fine.vol"));
sol1 = femBemCoupledSolve(model);
sol2 = femBemCoupledSolve(model, "MaterialCoef", 2);

[M, ~] = SurfaceP1Space(model.surface).mass();
area = sum(model.surface.areas());
% discretely the exterior part is invariant only up to the (small)
% variation of the computed trace: measured 1.5e-6 on this fixture.
verifyEqual(testCase, sum(M * sol2.trace) / area, ...
    sum(M * sol1.trace) / area, "AbsTol", 1e-4);
verifyEqual(testCase, sol2.totalExteriorFlux, sol1.totalExteriorFlux, ...
    "RelTol", 1e-4);
% interior deviation from the trace halves when c doubles
dev1 = max(sol1.u) - sum(M * sol1.trace) / area;
dev2 = max(sol2.u) - sum(M * sol2.trace) / area;
verifyEqual(testCase, dev2, dev1 / 2, "RelTol", 2e-2);
end


function surface = fixtureSurface(name)
mesh = VolMesh(fixturePath(name));
surface = mesh.boundary();
end


function path = fixturePath(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
path = string(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end
