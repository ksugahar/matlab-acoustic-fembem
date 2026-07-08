function tests = testCurvedPanelBem
%TESTCURVEDPANELBEM Fast checks for the isoparametric curved-panel BEM lane.
%
% Curved (quadratic-isoparametric) panels remove the O(h^2) straight-panel
% faceting error that caps the flat BEM (testHelmholtzScattering documents the
% "faceted-geometry dominated" deviation).  These are the quick regression
% gates; the heavy convergence A/B lives in
% validation_test/testCurvedPanelSphereConvergence.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testFlatProjectionReducesToSurfaceQuadrature(testCase)
% default Projection = @(X) X must reproduce SurfaceQuadrature bit-for-bit:
% the quadratic map with straight edge nodes degenerates to the affine one.
surface = fixtureSurface("unit_sphere_coarse.vol");
for order = [1 3 7]
    fq = CurvedPanelQuadrature(surface, order);
    sq = SurfaceQuadrature(surface, order);
    verifyLessThan(testCase, max(abs(fq.points(:) - sq.points(:))), 1e-12);
    verifyLessThan(testCase, max(abs(fq.weights - sq.weights)), 1e-12);
end
end


function testCurvedSurfaceAreaBeatsFlat(testCase)
% curving the sphere panels drives the total area to 4*pi*R^2 far faster.
surface = fixtureSurface("unit_sphere_coarse.vol");
proj = CurvedPanelQuadrature.sphereProjection(1.0);
flatArea = sum(CurvedPanelQuadrature(surface, 7).weights);
curvedArea = sum(CurvedPanelQuadrature(surface, 7, proj).weights);
areaTrue = 4 * pi;
verifyLessThan(testCase, abs(curvedArea - areaTrue) / areaTrue, 1e-3);
verifyLessThan(testCase, abs(curvedArea - areaTrue), 0.1 * abs(flatArea - areaTrue));
end


function testCurvedLaplaceCapacitanceBeatsFlat(testCase)
% potential 1 on the unit sphere -> total charge = capacitance = 4*pi*R.
surface = fixtureSurface("unit_sphere_coarse.vol");
proj = CurvedPanelQuadrature.sphereProjection(1.0);
g = ones(size(surface.vtx, 1), 1);
solF = curvedSingleLayerDirichletSolve(surface, g);
solC = curvedSingleLayerDirichletSolve(surface, g, "Projection", proj);
verifyEqual(testCase, solF.status, "ok");
verifyEqual(testCase, solC.status, "ok");
verifyTrue(testCase, isreal(solC.q));
eF = abs(solF.totalCharge - 4*pi) / (4*pi);
eC = abs(solC.totalCharge - 4*pi) / (4*pi);
verifyLessThan(testCase, eC, 1e-3);
verifyLessThan(testCase, eC, eF);
end


function testCurvedSoftSphereScatterBeatsFlat(testCase)
% sound-soft plane-wave scattering vs the analytic partial-wave series:
% the curved lane is well inside 1e-2 and clearly beats the flat lane.
surface = fixtureSurface("unit_sphere_coarse.vol");
proj = CurvedPanelQuadrature.sphereProjection(1.0);
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
k = 2.0;
g = -exp(1i * k * surface.vtx(:, 3));
solF = curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k);
solC = curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k, "Projection", proj);
verifyFalse(testCase, isreal(solC.q));
ref = softSphereScattering(k, 1.0, probes).scattered;
eF = max(abs(solF.potentialAt(probes) - ref) ./ abs(ref));
eC = max(abs(solC.potentialAt(probes) - ref) ./ abs(ref));
verifyLessThan(testCase, eC, 0.02);
verifyLessThan(testCase, eC, 0.5 * eF);
end


function surface = fixtureSurface(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
file = fullfile(repoRoot, "fixtures", "mesh_topology", name);
surface = VolMesh(file).boundary();
end
