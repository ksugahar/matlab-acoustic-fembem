function tests = testNgsolveBemCrossCheck
%TESTNGSOLVEBEMCROSSCHECK Standing regression against ngsolve.bem references.
%
% The committed validation/data/ngbem_reference_*.mat files hold dense V/K/M
% assembled by NGSolve's ngsolve.bem (Sauter-Schwab, intorder 16) on the
% sphere fixtures; regenerate with validation/exportNgsolveBemReference.py.
% Unlike the Gypsilab cross-check (which needs the external Gypsilab tree at
% run time) this one runs from the committed artifacts alone, so it lives in
% the test suite. Bands locked from the 2026-07-03 measurements.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "validation"));
end


function testCoarseSphereMatchesNgsolveBem(testCase)
report = crossCheck("unit_sphere_coarse");
% measured: mass 1.3e-16, V 3.80e-4, K 3.40e-3, capacitance 1.64e-4
verifyEqual(testCase, report.status, "ok");
verifyLessThan(testCase, report.massRelDiff, 1e-12);
verifyLessThan(testCase, report.operatorVRelDiff, 1e-3);
verifyLessThan(testCase, report.operatorKRelDiff, 8e-3);
verifyLessThan(testCase, report.capacitanceRelDiff, 1e-3);
end


function testFineSphereMatchesNgsolveBem(testCase)
report = crossCheck("unit_sphere_fine");
% measured: mass 1.2e-16, V 2.85e-4, K 3.24e-3, capacitance 4.28e-5;
% the fine fixture also pins that NGSolve's interior-vertex rows are zero
% (18 interior nodes dropped by the volNodeIds restriction).
verifyEqual(testCase, report.status, "ok");
verifyLessThan(testCase, report.massRelDiff, 1e-12);
verifyLessThan(testCase, report.operatorVRelDiff, 1e-3);
verifyLessThan(testCase, report.operatorKRelDiff, 8e-3);
verifyLessThan(testCase, report.capacitanceRelDiff, 1e-3);
end


function testReferenceArtifactsAreConverged(testCase)
% the .mat itself records its intorder 12-vs-16 self-convergence; a stale
% or hand-edited artifact fails here before any operator comparison.
for base = ["unit_sphere_coarse", "unit_sphere_fine"]
    report = crossCheck(base);
    verifyLessThan(testCase, report.referenceIntorderConvergenceV, 1e-6);
    verifyLessThan(testCase, report.referenceIntorderConvergenceK, 1e-6);
    verifyEqual(testCase, report.referenceIntorder, 16);
end
end


function report = crossCheck(baseName)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", baseName + ".vol");
report = verifyGalerkinAgainstNgsolve(volFile);
end
