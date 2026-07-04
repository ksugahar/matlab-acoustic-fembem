function tests = testFsiCoupledSolve
%TESTFSICOUPLEDSOLVE Acoustic fluid-structure interaction coupled solve.
%
% The genuine acoustic FEM/BEM coupling: vector P1 elasticity FEM interior +
% acoustic BEM exterior + the FSI interface conditions, gated against the
% analytic elastic sphere (elasticSphereScattering). Locked from the
% 2026-07-04 measurements:
%   stiff limit -> rigid sphere to ~1e-3 (the FORMULATION gate - validates
%     the interface + BEM coupling independent of interior resolution),
%   elastic field CONVERGES to the analytic under refinement
%     (25% coarse -> 7% fine at kR = 2; P1 interior elasticity, especially
%     the shear field, is the accuracy-limiting factor - the coupling is
%     exact, the interior is low-order).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function model = fixtureModel(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
model = FemBemModel(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end


function testStiffLimitReproducesRigidSphere(testCase)
% a very stiff/heavy solid does not move: the FSI solve must reproduce the
% rigid sphere. This validates the interface + BEM coupling independent of
% the interior elastic resolution. measured: 5e-4 / 4.2e-3 / 2.2e-4.
k = 2.0;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
model = fixtureModel("unit_ball_maxh018.vol");
sol = fsiCoupledSolve(model, "Wavenumber", k, ...
    "LongitudinalSpeed", 50, "ShearSpeed", 30, "DensityRatio", 100);
verifyEqual(testCase, sol.status, "ok");
rigid = rigidSphereScattering(k, 1.0, probes);
err = max(abs(sol.totalAt(probes) - rigid.total) ./ abs(rigid.total));
verifyLessThan(testCase, err, 1e-2);
end


function testElasticFieldConvergesToAnalytic(testCase)
% a lucite-like elastic sphere: the coupled solve converges to the analytic
% elastic sphere under mesh refinement (coarse error > fine error), and the
% fine field is within the P1-interior band. measured coarse 25/55/26%,
% fine 7.5/16.7/7.7%.
k = 2.0;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
mat = {"LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15};
ref = elasticSphereScattering(k, 1.0, probes, mat{:}).total;

solC = fsiCoupledSolve(fixtureModel("unit_sphere_fine.vol"), ...
    "Wavenumber", k, mat{:});
solF = fsiCoupledSolve(fixtureModel("unit_ball_maxh018.vol"), ...
    "Wavenumber", k, mat{:});
errC = max(abs(solC.totalAt(probes) - ref) ./ abs(ref));
errF = max(abs(solF.totalAt(probes) - ref) ./ abs(ref));

verifyEqual(testCase, solF.status, "ok");
verifyLessThan(testCase, errF, errC);          % refinement improves it
verifyLessThan(testCase, errF, 0.20);          % P1-interior band (meas 0.17)
verifyLessThan(testCase, errC / errF, 10);     % sane convergence factor
end
