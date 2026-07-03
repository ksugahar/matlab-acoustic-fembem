function tests = testGalerkinSingleLayer
%TESTGALERKINSINGLELAYER Galerkin BEM operator and the stage-4 sphere BVP.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testQuadratureRefinementConvergesAndSymmetrizes(testCase)
% the true Galerkin matrix is symmetric; the one-sided semi-analytic
% assembly must approach it as the test quadrature is refined.
surface = fixtureSurface("closed_tetra_surface_manifold.vol");
asym = zeros(1, 3);
orders = [1 3 7];
for k = 1:3
    op = GalerkinSingleLayer(surface, "QuadratureOrder", orders(k));
    asym(k) = norm(op.matrix - op.matrix.', "fro") / norm(op.matrix, "fro");
end
verifyTrue(testCase, asym(2) < asym(1) && asym(3) < asym(2));
verifyLessThan(testCase, asym(3), 2e-4);
end


function testLowFrequencyLimitMatchesLaplace(testCase)
surface = fixtureSurface("closed_tetra_surface_manifold.vol");
op0 = GalerkinSingleLayer(surface);
opk = GalerkinSingleLayer(surface, "Wavenumber", 1e-9);
relDiff = norm(opk.matrix - op0.matrix, "fro") / norm(op0.matrix, "fro");
verifyLessThan(testCase, relDiff, 1e-8);
verifyTrue(testCase, isreal(op0.matrix));
end


function testSymmetrizedOperatorIsPositiveDefinite(testCase)
% the Laplace single layer on a closed surface is positive definite.
surface = fixtureSurface("unit_sphere_coarse.vol");
op = GalerkinSingleLayer(surface);
ev = eig((op.matrix + op.matrix.') / 2);
verifyGreaterThan(testCase, min(ev), 0);
end


function testSphereCapacitanceMatchesAnalyticBand(testCase)
% stage-4 BVP: exterior Dirichlet g = 1 on the coarse unit sphere.
% Analytic capacitance is 4*pi; the faceted mesh sits ~2.9% low
% (geometry-dominated: gss 3 vs 7 moves it by 0.01%). Band locked from
% the 2026-07-03 measurement and the Gypsilab cross-check (rel diff 1.4e-5).
surface = fixtureSurface("unit_sphere_coarse.vol");
sol = singleLayerDirichletSolve(surface, ones(size(surface.vtx, 1), 1));

verifyEqual(testCase, sol.status, "ok");
verifyEqual(testCase, sol.kind, "laplace_single_layer_exterior_dirichlet_solve");
verifyLessThan(testCase, sol.residualNorm, 1e-12);
verifyEqual(testCase, sol.totalCharge, 12.2046, "AbsTol", 0.05);   % golden band
verifyLessThan(testCase, abs(sol.totalCharge - 4*pi) / (4*pi), 0.05);
end


function testFarFieldDecaysLikeAPointCharge(testCase)
surface = fixtureSurface("unit_sphere_coarse.vol");
sol = singleLayerDirichletSolve(surface, ones(size(surface.vtx, 1), 1));
points = [3 0 0; 0 0 5];
u = sol.potentialAt(points);
uPoint = sol.totalCharge ./ (4 * pi * [3; 5]);
verifyEqual(testCase, u, uPoint, "RelTol", 5e-3);
end


function testRejectsWrongBoundaryLength(testCase)
surface = fixtureSurface("unit_sphere_coarse.vol");
verifyError(testCase, @() singleLayerDirichletSolve(surface, [1; 2]), ...
    "singleLayerDirichletSolve:boundary");
end


function testSurfaceQuadratureIntegratesExactly(testCase)
% the rules must integrate constants (all orders) and P1 (order >= 3)
% exactly: sum of weights = area, integral of each basis = area/3 per tri.
surface = fixtureSurface("closed_tetra_surface_manifold.vol");
totalArea = sum(surface.areas());
for order = [1 3 7]
    quad = SurfaceQuadrature(surface, order);
    verifyEqual(testCase, sum(quad.weights), totalArea, "RelTol", 1e-13);
    basisIntegrals = full(sum(quad.weightedBasis(), 1)).';
    [M, ~] = SurfaceP1Space(surface).mass();
    verifyEqual(testCase, basisIntegrals, full(sum(M, 2)), "RelTol", 1e-12);
end
end


function surface = fixtureSurface(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
mesh = VolMesh(string(fullfile(repoRoot, "fixtures", "mesh_topology", name)));
surface = mesh.boundary();
end
