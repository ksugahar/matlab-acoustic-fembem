function tests = testCoulombGauge
%testCoulombGauge Tests for MQS Coulomb-gauge postprocess gates.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testMqCoulombGaugePackagePasses(testCase)
artifacts = baseArtifacts();

result = mqCoulombGaugePostprocessPackage(artifacts, ...
    "ExpectedCaseId", "coil_mqs_001", ...
    "ExpectedMeshId", "unit_tet_mesh", ...
    "ExpectedFrequencyHz", 1.0e6);

verifyEqual(testCase, result.kind, "mqs_coulomb_gauge_efield_postprocess_package");
verifyEqual(testCase, result.policy, "readable_mqs_coulomb_gauge_efield_postprocess_gate");
verifyEqual(testCase, result.status, "ok");
verifyTrue(testCase, result.checks.requiredKindsPresent);
verifyTrue(testCase, result.checks.mqsAPhiFormulationRecorded);
verifyTrue(testCase, result.checks.coulombGaugeConditionRecorded);
verifyTrue(testCase, result.checks.spatialPotentialBcRecorded);
verifyTrue(testCase, result.checks.electricFieldUnitRecorded);
verifyTrue(testCase, result.checks.validityEnvelopeOk);
verifyEqual(testCase, result.caseIds, "coil_mqs_001");
verifyEqual(testCase, result.meshIds, "unit_tet_mesh");
verifyEqual(testCase, result.frequenciesHz, 1.0e6);
end


function testRejectsStaleFrequency(testCase)
artifacts = baseArtifacts();
artifacts(4).frequencyHz = 1.2e6;

result = mqCoulombGaugePostprocessPackage(artifacts);

verifyEqual(testCase, result.status, "needs_attention");
verifyFalse(testCase, result.checks.frequenciesMatch);
end


function testRejectsMissingBoundaryConditionSource(testCase)
artifacts = baseArtifacts();
artifacts(3).boundaryConditionSource = "";

result = mqCoulombGaugePostprocessPackage(artifacts);

verifyEqual(testCase, result.status, "needs_attention");
verifyFalse(testCase, result.checks.spatialPotentialBcRecorded);
end


function testRejectsInvalidValidityEnvelope(testCase)
artifacts = baseArtifacts();
artifacts(5).frequencyRatioToFullwaveLimit = 0.5;

result = mqCoulombGaugePostprocessPackage(artifacts, ...
    "MaxFrequencyRatioToFullwave", 0.1);

verifyEqual(testCase, result.status, "needs_attention");
verifyFalse(testCase, result.checks.validityEnvelopeOk);
end


function artifacts = baseArtifacts()
template = struct( ...
    "kind", "", ...
    "caseId", "coil_mqs_001", ...
    "meshId", "unit_tet_mesh", ...
    "frequencyHz", 1.0e6, ...
    "path", "", ...
    "status", "ok", ...
    "gatePolicy", "", ...
    "formulation", "", ...
    "gaugeCondition", "", ...
    "boundaryConditionSource", "", ...
    "EUnit", "", ...
    "frequencyRatioToFullwaveLimit", NaN, ...
    "dominantInductive", false, ...
    "comparisonReference", "");

artifacts = repmat(template, 1, 5);

artifacts(1).kind = "mqs_solution";
artifacts(1).path = "slot160_mqs_solution.json";
artifacts(1).gatePolicy = "mqs_a_phi_solution";
artifacts(1).formulation = "A_phi";

artifacts(2).kind = "coulomb_gauge";
artifacts(2).path = "slot160_coulomb_gauge.json";
artifacts(2).gatePolicy = "coulomb_gauge_postprocess";
artifacts(2).gaugeCondition = "div_A_zero";

artifacts(3).kind = "spatial_potential";
artifacts(3).path = "slot160_spatial_potential.json";
artifacts(3).gatePolicy = "electrostatic_potential_postprocess";
artifacts(3).boundaryConditionSource = "conductor_surface_potential_from_coulomb_gauge";

artifacts(4).kind = "electric_field";
artifacts(4).path = "slot160_efield.json";
artifacts(4).gatePolicy = "efield_gradient_recovery";
artifacts(4).EUnit = "V_per_m";

artifacts(5).kind = "validity_envelope";
artifacts(5).path = "slot160_validity.json";
artifacts(5).gatePolicy = "mqs_darwin_fullwave_validity_envelope";
artifacts(5).frequencyRatioToFullwaveLimit = 0.02;
artifacts(5).dominantInductive = true;
artifacts(5).comparisonReference = "Darwin_full_wave";
end
