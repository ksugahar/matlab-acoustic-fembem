function tests = testAdjointScalingManifest
tests = functiontests(localfunctions);
end


function rows = sampleRows()
counts = [4 8 16];
errors = [2e-10 5e-10 4e-9];
ratios = [2.1 4.5 10.8];
rows = repmat(struct(), numel(counts), 1);
for k = 1:numel(counts)
    rows(k).designVariableCount = counts(k);
    rows(k).adjointSolves = 1;
    rows(k).gradientCheckRelativeError = errors(k);
    rows(k).forwardAffineResidual = 1e-17;
    rows(k).plusAscentObjectiveRatio = 1.01;
    rows(k).minusAscentObjectiveRatio = 0.99;
    rows(k).fiftyStepObjectiveRatio = ratios(k);
    rows(k).fiftyStepMonotone = true;
end
end


function testAcceptsOneSolveAndSignedAscent(testCase)
report = adjointScalingManifest(sampleRows());
verifyTrue(testCase, report.ok);
verifyEqual(testCase, report.adjointSolves, [1 1 1]);
out = evalc("acoustic_fembem.check_adjoint_scaling(jsonencode(sampleRows()), 1e-6, 1e-10, 1.0)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_adjoint_scaling");
end


function testRejectsVariableCostAndWrongDirection(testCase)
rows = sampleRows();
rows(3).adjointSolves = 16;
rows(2).plusAscentObjectiveRatio = 0.98;
report = adjointScalingManifest(rows);
verifyFalse(testCase, report.ok);
verifyFalse(testCase, report.checks.oneAdjointSolveForEverySize);
verifyFalse(testCase, report.checks.positiveDirectionRaisesObjective);
end
