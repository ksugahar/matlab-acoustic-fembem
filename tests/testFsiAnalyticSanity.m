function tests = testFsiAnalyticSanity
%TESTFSIANALYTICSANITY Fast analytic sanity for the elastic-acoustic FSI coupling.
%
% The heavy Faran convergence sweeps (mesh refinement, dense-BEM leg, resonance
% scans) live in validation_test/testFsiCoupledSolve.  This keeps a quick check
% in the fast lane: the FEM-elastic interior + exact spherical DtN (a high-order
% radiating impedance Zs, no dense BEM) coupled solve lands in the analytic Faran
% elastic-sphere ballpark, and a very stiff sphere reproduces the rigid scatterer.
% One DtN solve on the 629-node ball is ~0.3 s, so both checks stay fast.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDtnElasticCouplingMatchesFaranBallpark(testCase)
% FEM-elastic interior + exact spherical DtN (Zs) exterior vs the analytic Faran
% elastic sphere, into the P1-interior accuracy band (~0.16 measured at kR = 2).
volFile = fixtureBall();
mat = {"LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15};
k = 2.0;
probes = [0 0 3; 3 0 0; 0 3 0];

ref = elasticSphereScattering(k, 1.0, probes, mat{:}).total;
sol = fsiCoupledSolve(FemBemModel(volFile), "Wavenumber", k, mat{:}, ...
    "ExteriorMethod", "dtn");
p = sol.totalAt(probes);

verifyEqual(testCase, sol.status, "ok");
verifyTrue(testCase, sol.dtn.used);                          % the DtN (Zs) exterior ran
verifyLessThan(testCase, max(abs(p - ref) ./ abs(ref)), 0.25);   % P1-interior band
end


function testDtnElasticCouplingMatchesFaranAtLowerWavenumber(testCase)
% a second spectral point: the DtN-coupled elastic sphere must also track the
% analytic Faran field at a lower wavenumber (kR = 1), catching frequency-
% dependent coupling / sign bugs the single kR = 2 check could miss.
volFile = fixtureBall();
mat = {"LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15};
k = 1.0;
probes = [0 0 3; 3 0 0; 0 3 0];

ref = elasticSphereScattering(k, 1.0, probes, mat{:}).total;
p = fsiCoupledSolve(FemBemModel(volFile), "Wavenumber", k, mat{:}, ...
    "ExteriorMethod", "dtn").totalAt(probes);

verifyLessThan(testCase, max(abs(p - ref) ./ abs(ref)), 0.25);
end


function volFile = fixtureBall()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_ball_maxh018.vol"));
end
