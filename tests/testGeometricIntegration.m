function tests = testGeometricIntegration
%testGeometricIntegration Tests for geometric integration teaching.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testHarmonicOscillatorEnergyReport(testCase)
report = geometricIntegratorEnergyReport(0.02, 1000, 1.0);

verifyEqual(testCase, report.kind, "geometric_integrator_energy_report");
verifyEqual(testCase, report.policy, "readable_geometric_time_integration_energy_gate");
verifyTrue(testCase, report.pass);
verifyEqual(testCase, numel(report.method_rows), 3);
verifyEqual(testCase, report.steps, 1000);
verifyEqual(testCase, report.step_size_s, 0.02, "AbsTol", 0);
verifyEqual(testCase, report.omega_rad_per_s, 1.0, "AbsTol", 0);

verifyGreaterThan(testCase, report.explicit_euler.max_rel_energy_drift, 0.40);
verifyLessThan(testCase, report.symplectic_euler.max_rel_energy_drift, 0.02);
verifyLessThan(testCase, report.implicit_midpoint.max_rel_energy_drift, 1e-10);
verifyGreaterThan(testCase, report.explicit_to_geometric_drift_ratio, 20);
verifyTrue(testCase, report.checks.symplecticEulerBounded);
verifyTrue(testCase, report.checks.implicitMidpointPreservesQuadraticEnergy);
verifyTrue(testCase, report.checks.explicitEulerIsNegativeControl);
end


function testRejectsBadInputs(testCase)
verifyError(testCase, @() geometricIntegratorEnergyReport(0, 100, 1), ...
    "geometricIntegratorEnergyReport:stepSize");
verifyError(testCase, @() geometricIntegratorEnergyReport(0.02, 10.5, 1), ...
    "geometricIntegratorEnergyReport:steps");
verifyError(testCase, @() geometricIntegratorEnergyReport(0.02, 100, 0), ...
    "geometricIntegratorEnergyReport:omega");
end
