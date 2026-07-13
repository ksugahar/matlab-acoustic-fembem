function tests = testFsiScatteringInvariantsGate
tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);
end


function testAcceptsMeasuredLosslessFsiEvidence(testCase)
report = acoustic_fembem.fsi_scattering_invariants_gate( ...
    6.370374314448229e-4, 8.140261530724266e-4, 1.3152264144524783e-2, ...
    8.596141376286232e-15, true, "exp(-i omega t)");
verifyTrue(testCase, report.ok);
verifyTrue(testCase, report.checks.opticalTheoremEnergyClosure);
verifyTrue(testCase, report.checks.p1BemHighOrderDtnAgreement);
end


function testRejectsEnergyAndExteriorRegression(testCase)
report = acoustic_fembem.fsi_scattering_invariants_gate( ...
    1e-3, 0.08, 0.07, 1e-14, true, "exp(-i omega t)");
verifyFalse(testCase, report.ok);
verifyFalse(testCase, report.checks.opticalTheoremEnergyClosure);
verifyFalse(testCase, report.checks.p1BemHighOrderDtnAgreement);
end


function testMcpWrapperRejectsAmbiguousPhysicalClaim(testCase)
out = evalc("acoustic_fembem.check_fsi_scattering_invariants(1e-3, 1e-3, 1e-2, 1e-14, false, ""unspecified"")");
decoded = jsondecode(out);
verifyFalse(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_fsi_scattering_invariants");
verifyFalse(testCase, decoded.result.checks.losslessMaterialDeclared);
verifyFalse(testCase, decoded.result.checks.outgoingTimeConventionExplicit);
end
