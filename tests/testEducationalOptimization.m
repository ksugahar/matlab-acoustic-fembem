function tests = testEducationalOptimization
%TESTEDUCATIONALOPTIMIZATION Tests for readable optimization gates.

tests = functiontests(localfunctions);
end


function testExactLeastSquaresAndGradientCheck(testCase)
A = [ ...
    1 0
    0 1
    1 1
    2 -1];
xTrue = [2; -1];
b = A * xTrue;
x0 = [0.25; -0.5];

result = educationalQuadraticLeastSquares(A, b, "Initial", x0);

verifyEqual(testCase, result.kind, "educational_quadratic_least_squares");
verifyEqual(testCase, result.policy, "readable_matlab_optimization_gate_not_optuna_owned");
verifyEqual(testCase, result.x, xTrue, "AbsTol", 1e-12);
verifyLessThan(testCase, result.residualNorm, 1e-12);
verifyLessThan(testCase, result.normalEquationResidual, 1e-12);
verifyTrue(testCase, result.gradientCheck.passed);
verifyLessThan(testCase, result.gradientCheck.maxAbsError, 1e-8);
verifyGreaterThan(testCase, result.objectiveAtInitial, result.objective);
end


function testOverdeterminedLeastSquaresMatchesBackslash(testCase)
A = [ ...
    1 0
    1 1
    1 2
    1 3
    1 4];
b = [1.0; 1.8; 3.2; 3.9; 5.1];
x0 = [0; 0];

result = educationalQuadraticLeastSquares(A, b, ...
    "Initial", x0, ...
    "FiniteDifferenceStep", 1e-6, ...
    "GradientTolerance", 1e-7);

verifyEqual(testCase, result.x, A \ b, "AbsTol", 1e-12);
verifyLessThan(testCase, result.normalEquationResidual, 1e-12);
verifyTrue(testCase, result.gradientCheck.passed);
verifyGreaterThan(testCase, result.objectiveAtInitial, result.objective);
verifyEqual(testCase, result.matrix.rows, 5);
verifyEqual(testCase, result.matrix.cols, 2);
verifyEqual(testCase, result.matrix.rank, 2);
end


function testRejectsWrongInitialLength(testCase)
A = eye(2);
b = [1; 2];

verifyError(testCase, @() educationalQuadraticLeastSquares(A, b, "Initial", 0), ...
    "educationalQuadraticLeastSquares:initial");
end


function testBoxConstrainedLeastSquaresClampsUnconstrainedSolution(testCase)
A = eye(2);
b = [2; -1];
lower = [0; 0];
upper = [1; 3];

result = educationalBoxConstrainedLeastSquares(A, b, lower, upper, ...
    "Initial", [0; 0], ...
    "StepSize", 1.0, ...
    "MaxIterations", 5);

verifyEqual(testCase, result.kind, "educational_box_constrained_least_squares");
verifyEqual(testCase, result.policy, "readable_box_projected_gradient_gate_not_optuna_owned");
verifyEqual(testCase, result.x, [1; 0], "AbsTol", 1e-12);
verifyEqual(testCase, result.gradient, [-1; 1], "AbsTol", 1e-12);
verifyTrue(testCase, result.activeUpper(1));
verifyFalse(testCase, result.activeUpper(2));
verifyFalse(testCase, result.activeLower(1));
verifyTrue(testCase, result.activeLower(2));
verifyLessThan(testCase, result.maxKktResidual, 1e-12);
verifyLessThan(testCase, result.projectedGradientResidual, 1e-12);
verifyGreaterThan(testCase, result.objectiveHistory(1), result.objective);
verifyTrue(testCase, result.objectiveMonotone);
verifyTrue(testCase, result.gradientCheck.passed);
end


function testBoxConstrainedLeastSquaresInteriorSolution(testCase)
A = eye(2);
b = [0.25; 0.75];
lower = [0; 0];
upper = [1; 1];

result = educationalBoxConstrainedLeastSquares(A, b, lower, upper, ...
    "Initial", [0; 0], ...
    "StepSize", 1.0);

verifyEqual(testCase, result.x, b, "AbsTol", 1e-12);
verifyFalse(testCase, any(result.activeLower));
verifyFalse(testCase, any(result.activeUpper));
verifyLessThan(testCase, result.maxKktResidual, 1e-12);
verifyTrue(testCase, result.gradientCheck.passed);
end


function testBoxConstrainedLeastSquaresRejectsInvertedBounds(testCase)
A = eye(1);
b = 1;

verifyError(testCase, @() educationalBoxConstrainedLeastSquares(A, b, 2, 1), ...
    "educationalBoxConstrainedLeastSquares:bounds");
end


function testFemBemTraceLeastSquaresFitsBoundaryTrace(testCase)
path = writeFixture(testCase, tetVolText());
model = volFemBemModel(path);
target = [1; 2; 3; 4];

result = educationalFemBemTraceLeastSquares(model, target);

verifyEqual(testCase, result.kind, "educational_fem_bem_trace_least_squares");
verifyEqual(testCase, result.policy, "readable_trace_optimization_gate_first_order_tri_tet");
verifyEqual(testCase, result.trace, target, "AbsTol", 1e-12);
verifyLessThan(testCase, result.traceResidualNorm, 1e-12);
verifyLessThan(testCase, result.normalEquationResidual, 1e-12);
verifyTrue(testCase, result.gradientCheck.passed);
verifyGreaterThan(testCase, result.objectiveAtInitial, result.objective);
verifyEqual(testCase, result.matrix.traceRows, 4);
verifyEqual(testCase, result.matrix.femUnknowns, 4);
end


function testFemBemTraceLeastSquaresLeavesInteriorNodeUnforced(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());
model = volFemBemModel(path);
target = [10; 20; 30; 40];

result = educationalFemBemTraceLeastSquares(model, target, ...
    "Tikhonov", 1e-3, ...
    "Initial", ones(5, 1), ...
    "GradientTolerance", 1e-6);

verifyEqual(testCase, numel(result.u), 5);
verifyLessThan(testCase, abs(result.u(5)), 1e-12);
verifyTrue(testCase, result.gradientCheck.passed);
verifyGreaterThan(testCase, result.objectiveAtInitial, result.objective);
verifyEqual(testCase, result.matrix.traceRows, 4);
verifyEqual(testCase, result.matrix.femUnknowns, 5);
end


function testFemBemTikhonovPathShowsResidualNormTradeoff(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());
model = volFemBemModel(path);
target = [10; 20; 30; 40];

regularization = educationalFemBemTikhonovPath(model, target, [0; 1e-3; 1e-1; 1], ...
    "Initial", ones(5, 1), ...
    "GradientTolerance", 1e-6);

verifyEqual(testCase, regularization.kind, "educational_fem_bem_tikhonov_path");
verifyEqual(testCase, regularization.policy, "readable_trace_regularization_path_first_order_tri_tet");
verifyEqual(testCase, regularization.status, "ok");
verifyTrue(testCase, regularization.checks.alphasSorted);
verifyTrue(testCase, regularization.checks.solutionNormNonincreasing);
verifyTrue(testCase, regularization.checks.traceResidualNondecreasing);
verifyTrue(testCase, regularization.checks.weightedResidualNondecreasing);
verifyTrue(testCase, regularization.checks.allGradientChecksPassed);
verifyEqual(testCase, regularization.rows(1).solver, "minimum_norm_pinv_rank_deficient");
verifyLessThan(testCase, regularization.solutionNorms(end), regularization.solutionNorms(1));
verifyGreaterThan(testCase, regularization.traceResidualNorms(end), regularization.traceResidualNorms(1));
verifyLessThan(testCase, regularization.rows(end).gradientCheckMaxAbsError, 1e-6);
end


function testFemBemLcurveCornerSelectsInteriorAlpha(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());
model = volFemBemModel(path);
target = [10; 20; 30; 40];

corner = educationalFemBemLcurveCorner(model, target, [1e-4; 1e-3; 1e-2; 1e-1; 1], ...
    "Initial", ones(5, 1), ...
    "GradientTolerance", 1e-6);

verifyEqual(testCase, corner.kind, "educational_fem_bem_lcurve_corner");
verifyEqual(testCase, corner.policy, "readable_lcurve_corner_from_trace_regularization_path");
verifyEqual(testCase, corner.status, "ok");
verifyTrue(testCase, corner.checks.pathOk);
verifyTrue(testCase, corner.checks.interiorSelected);
verifyTrue(testCase, corner.checks.finiteCurvature);
verifyTrue(testCase, corner.checks.positiveCornerCurvature);
verifyGreaterThan(testCase, corner.selectedIndex, 1);
verifyLessThan(testCase, corner.selectedIndex, numel(corner.alphas));
verifyEqual(testCase, corner.selectedAlpha, corner.alphas(corner.selectedIndex));
verifyEqual(testCase, corner.selectedTraceResidualNorm, ...
    corner.path.traceResidualNorms(corner.selectedIndex), "AbsTol", 1e-14);
verifyEqual(testCase, corner.selectedSolutionNorm, ...
    corner.path.solutionNorms(corner.selectedIndex), "AbsTol", 1e-14);
verifyEqual(testCase, corner.curvature(1), 0);
verifyEqual(testCase, corner.curvature(end), 0);
end


function testFemBemMorozovDiscrepancySelectsNoiseMatchedAlpha(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());
model = volFemBemModel(path);
target = [10; 20; 30; 40];
alphas = [0; 1e-3; 1e-2; 1e-1; 1];

regularization = educationalFemBemTikhonovPath(model, target, alphas, ...
    "Initial", ones(5, 1), ...
    "GradientTolerance", 1e-6);
noiseNorm = regularization.weightedTraceResiduals(4);

choice = educationalFemBemMorozovDiscrepancy(model, target, alphas, noiseNorm, ...
    "Initial", ones(5, 1), ...
    "GradientTolerance", 1e-6);

verifyEqual(testCase, choice.kind, "educational_fem_bem_morozov_discrepancy");
verifyEqual(testCase, choice.policy, "readable_morozov_discrepancy_from_trace_regularization_path");
verifyEqual(testCase, choice.status, "ok");
verifyTrue(testCase, choice.checks.pathOk);
verifyTrue(testCase, choice.checks.residualsNondecreasing);
verifyTrue(testCase, choice.checks.noiseBracketed);
verifyTrue(testCase, choice.checks.selectedMinimizesDiscrepancy);
verifyEqual(testCase, choice.selectedIndex, 4);
verifyEqual(testCase, choice.selectedAlpha, alphas(4), "AbsTol", 1e-14);
verifyEqual(testCase, choice.selectedResidualNorm, noiseNorm, "AbsTol", 1e-12);
verifyEqual(testCase, choice.lowerBracketIndex, 4);
verifyEqual(testCase, choice.upperBracketIndex, 4);
end


function testTouchstonePortMetadataFreezesOptionLineBeforeRows(testCase)
result = educationalTouchstonePortMetadata(["P1" "P2"], ...
    "RequiredPorts", ["P1" "P2"], ...
    "PortOrder", ["P1" "P2"], ...
    "NetworkParameter", "S", ...
    "DataFormat", "MA", ...
    "ExpectedDataFormat", "MA", ...
    "FrequencyUnit", "GHz", ...
    "ExpectedFrequencyUnit", "GHz", ...
    "Z0", 50, ...
    "ExpectedZ0", 50);

verifyEqual(testCase, result.kind, "educational_touchstone_port_metadata");
verifyEqual(testCase, result.policy, "readable_touchstone_metadata_before_row_values");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.ports, ["P1"; "P2"]);
verifyTrue(testCase, result.checks.portOrderMatchesExpected);
verifyTrue(testCase, result.checks.touchstoneFormatMatchesExpected);
verifyTrue(testCase, result.checks.referenceImpedanceMatchesExpected);

swapped = educationalTouchstonePortMetadata(["P2" "P1"], ...
    "PortOrder", ["P1" "P2"], ...
    "NetworkParameter", "S", ...
    "DataFormat", "MA", ...
    "FrequencyUnit", "GHz", ...
    "Z0", 50);
verifyEqual(testCase, swapped.status, "needs_attention");
verifyFalse(testCase, swapped.checks.portOrderMatchesExpected);

missingZ0 = educationalTouchstonePortMetadata(["P1" "P2"], ...
    "PortOrder", ["P1" "P2"], ...
    "NetworkParameter", "S", ...
    "DataFormat", "MA", ...
    "FrequencyUnit", "GHz");
verifyEqual(testCase, missingZ0.status, "needs_attention");
verifyFalse(testCase, missingZ0.checks.referenceImpedanceRecorded);

wrongFormat = educationalTouchstonePortMetadata(["P1" "P2"], ...
    "PortOrder", ["P1" "P2"], ...
    "NetworkParameter", "S", ...
    "DataFormat", "DB", ...
    "ExpectedDataFormat", "MA", ...
    "FrequencyUnit", "GHz", ...
    "Z0", 50);
verifyEqual(testCase, wrongFormat.status, "needs_attention");
verifyFalse(testCase, wrongFormat.checks.touchstoneFormatMatchesExpected);
end


function testFarfieldPatternMetadataFreezesCutsUnitsAndPolarization(testCase)
result = educationalFarfieldPatternMetadata([0 90 180], [0 90], ...
    "FrequencyHz", 2.45e9, ...
    "ExpectedFrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RequiredPhiValuesDeg", [0 90], ...
    "RowCount", 6);

verifyEqual(testCase, result.kind, "educational_farfield_pattern_metadata");
verifyEqual(testCase, result.policy, "readable_farfield_metadata_before_pattern_values");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.thetaSpanDeg, 180, "AbsTol", 1e-12);
verifyEqual(testCase, result.expectedGridRows, 6);
verifyTrue(testCase, result.checks.requiredComponentsPresent);
verifyTrue(testCase, result.checks.requiredPhiValuesPresent);
verifyTrue(testCase, result.checks.normalizationMatchesExpected);

wrongBasis = educationalFarfieldPatternMetadata([0 90 180], [0 90], ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "linear", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RowCount", 6);
verifyEqual(testCase, wrongBasis.status, "needs_attention");
verifyFalse(testCase, wrongBasis.checks.polarizationBasisMatchesExpected);

missingNormalization = educationalFarfieldPatternMetadata([0 90 180], [0 90], ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RowCount", 6);
verifyEqual(testCase, missingNormalization.status, "needs_attention");
verifyFalse(testCase, missingNormalization.checks.normalizationRecorded);

narrowTheta = educationalFarfieldPatternMetadata([0 60 120], [0 90], ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RowCount", 6);
verifyEqual(testCase, narrowTheta.status, "needs_attention");
verifyFalse(testCase, narrowTheta.checks.thetaSpanOk);
end


function testFarfieldGainDirectivityEfficiencyChecksAcceptedPowerRows(testCase)
directivityDbi = 7.0;
eta = 0.65;
gainDbi = 10 * log10((10^(directivityDbi / 10)) * eta);

result = educationalFarfieldGainDirectivityEfficiency(gainDbi, directivityDbi, ...
    "RadiatedPowerW", 6.5, ...
    "AcceptedPowerW", 10.0, ...
    "Normalization", "accepted_power");

verifyEqual(testCase, result.kind, "educational_farfield_gain_directivity_efficiency");
verifyEqual(testCase, result.policy, "readable_farfield_gain_after_accepted_power_metadata");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.radiationEfficiency, eta, "AbsTol", 1e-14);
verifyEqual(testCase, result.efficiencyFromGainOverDirectivity, eta, "AbsTol", 1e-14);
verifyLessThan(testCase, result.gainRelativeError, 1e-12);
verifyTrue(testCase, result.checks.gainMatchesDirectivityTimesEfficiency);

tooHighGain = educationalFarfieldGainDirectivityEfficiency(directivityDbi + 0.2, directivityDbi, ...
    "RadiatedPowerW", 6.5, ...
    "AcceptedPowerW", 10.0, ...
    "Normalization", "accepted_power");
verifyEqual(testCase, tooHighGain.status, "needs_attention");
verifyFalse(testCase, tooHighGain.checks.gainNotAboveDirectivity);

missingNormalization = educationalFarfieldGainDirectivityEfficiency(gainDbi, directivityDbi, ...
    "RadiatedPowerW", 6.5, ...
    "AcceptedPowerW", 10.0);
verifyEqual(testCase, missingNormalization.status, "needs_attention");
verifyFalse(testCase, missingNormalization.checks.normalizationRecorded);

inconsistentPower = educationalFarfieldGainDirectivityEfficiency(gainDbi, directivityDbi, ...
    "RadiatedPowerW", 5.0, ...
    "AcceptedPowerW", 10.0, ...
    "Normalization", "accepted_power");
verifyEqual(testCase, inconsistentPower.status, "needs_attention");
verifyFalse(testCase, inconsistentPower.checks.gainMatchesDirectivityTimesEfficiency);
end


function testFarfieldLobeNotebookHandoffKeepsIdentityAndGainGate(testCase)
directivityDbi = 7.0;
eta = 0.65;
gainDbi = 10 * log10((10^(directivityDbi / 10)) * eta);
lobe = struct( ...
    "lobeId", "main", ...
    "frequencyHz", 2.45e9, ...
    "thetaDeg", 90.0, ...
    "phiDeg", 0.0, ...
    "polarizationBasis", "theta_phi", ...
    "normalization", "accepted_power", ...
    "gainUnit", "dBi", ...
    "directivityUnit", "dBi", ...
    "gainDbi", gainDbi, ...
    "directivityDbi", directivityDbi, ...
    "radiatedPowerW", 6.5, ...
    "acceptedPowerW", 10.0);

result = educationalFarfieldLobeNotebookHandoff([0 90 180], [0 90], lobe, ...
    "FrequencyHz", 2.45e9, ...
    "ExpectedFrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RequiredPhiValuesDeg", [0 90], ...
    "RowCount", 6);

verifyEqual(testCase, result.kind, "educational_farfield_lobe_notebook_handoff");
verifyEqual(testCase, result.policy, "readable_farfield_lobe_row_before_notebook_panel");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.lobeId, "main");
verifyEqual(testCase, result.gainDirectivity.radiationEfficiency, eta, "AbsTol", 1e-14);
verifyLessThan(testCase, result.gainDirectivity.gainRelativeError, 1e-12);
verifyTrue(testCase, result.checks.metadataOk);
verifyTrue(testCase, result.checks.gainDirectivityOk);
verifyTrue(testCase, result.checks.phiOnExportGrid);

missingLobe = rmfield(lobe, "lobeId");
missing = educationalFarfieldLobeNotebookHandoff([0 90 180], [0 90], missingLobe, ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RowCount", 6);
verifyEqual(testCase, missing.status, "needs_attention");
verifyFalse(testCase, missing.checks.lobeIdRecorded);

wrongPhi = lobe;
wrongPhi.phiDeg = 45.0;
wrongCut = educationalFarfieldLobeNotebookHandoff([0 90 180], [0 90], wrongPhi, ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RowCount", 6);
verifyEqual(testCase, wrongCut.status, "needs_attention");
verifyFalse(testCase, wrongCut.checks.phiOnExportGrid);

tooHigh = lobe;
tooHigh.gainDbi = directivityDbi + 0.2;
badGain = educationalFarfieldLobeNotebookHandoff([0 90 180], [0 90], tooHigh, ...
    "FrequencyHz", 2.45e9, ...
    "AngleUnit", "deg", ...
    "CoordinateSystem", "spherical", ...
    "PolarizationBasis", "theta_phi", ...
    "Quantity", "gain", ...
    "QuantityUnit", "dBi", ...
    "Normalization", "accepted_power", ...
    "FieldComponents", ["Etheta" "Ephi"], ...
    "RowCount", 6);
verifyEqual(testCase, badGain.status, "needs_attention");
verifyFalse(testCase, badGain.checks.gainDirectivityOk);
end


function testTouchstoneEquivalentCircuitKeepsReferenceImpedanceVisible(testCase)
result = educationalTouchstoneEquivalentCircuit(0, 0, ...
    "S12", 0, ...
    "S22", 0, ...
    "Z0", 50, ...
    "ComparisonZ0", 75);

verifyEqual(testCase, result.kind, "educational_touchstone_equivalent_circuit");
verifyEqual(testCase, result.policy, "readable_touchstone_z0_equivalent_circuit_gate");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.pi.yShunt1, 1/50, "AbsTol", 1e-14);
verifyEqual(testCase, result.pi.yShunt2, 1/50, "AbsTol", 1e-14);
verifyEqual(testCase, result.pi.ySeries, 0, "AbsTol", 1e-14);
verifyEqual(testCase, result.z.z11, 50, "AbsTol", 1e-14);
verifyEqual(testCase, result.t.zSeries1, 50, "AbsTol", 1e-14);
verifyEqual(testCase, result.t.zSeries2, 50, "AbsTol", 1e-14);
verifyEqual(testCase, result.t.zShunt, 0, "AbsTol", 1e-14);
verifyEqual(testCase, result.comparison.pi.yShunt1, 1/75, "AbsTol", 1e-14);
verifyEqual(testCase, result.comparison.z.z11, 75, "AbsTol", 1e-14);
verifyTrue(testCase, result.checks.equivalentValuesDependOnZ0);

active = educationalTouchstoneEquivalentCircuit(0, 1.1, "S12", 1.1, "S22", 0);
verifyEqual(testCase, active.status, "needs_attention");
verifyFalse(testCase, active.checks.sparameterPassivityOk);

verifyError(testCase, @() educationalTouchstoneEquivalentCircuit(0, 0, "ComparisonZ0", -75), ...
    "educationalTouchstoneEquivalentCircuit:comparisonZ0");
end


function testTouchstonePortMatchKeepsS11SeparateFromInsertionLoss(testCase)
result = educationalTouchstonePortMatch(1/3, ...
    "ReturnLossMinDb", 9.0, ...
    "VswrMax", 2.1);

verifyEqual(testCase, result.kind, "educational_touchstone_port_match");
verifyEqual(testCase, result.policy, "readable_one_port_match_quality_gate");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.reflectionCoefficient, 1/3, "AbsTol", 1e-14);
verifyEqual(testCase, result.vswr, 2.0, "AbsTol", 1e-14);
verifyEqual(testCase, result.returnLossDb, 9.542425094393248, "AbsTol", 1e-12);
verifyEqual(testCase, result.mismatchLossDb, -10 * log10(8/9), "AbsTol", 1e-14);
verifyEqual(testCase, result.transmittedPowerFraction, 8/9, "AbsTol", 1e-14);
verifyEqual(testCase, result.reflectedPowerFraction, 1/9, "AbsTol", 1e-14);
verifyTrue(testCase, result.checks.returnLossLimitOk);
verifyTrue(testCase, result.checks.vswrLimitOk);

matched = educationalTouchstonePortMatch(0);
verifyEqual(testCase, matched.vswr, 1.0, "AbsTol", 1e-14);
verifyTrue(testCase, isinf(matched.returnLossDb));
verifyEqual(testCase, matched.mismatchLossDb, 0.0, "AbsTol", 1e-14);

active = educationalTouchstonePortMatch(1.02);
verifyEqual(testCase, active.status, "needs_attention");
verifyFalse(testCase, active.checks.passiveReflectionOk);
verifyLessThan(testCase, active.transmittedPowerFraction, 0);
end


function testTouchstoneDesignFrequencyGridRequiresBracketDensity(testCase)
result = educationalTouchstoneDesignFrequencyGrid([0.95 0.99 1.01 1.05], 1.0, ...
    "FrequencyUnit", "GHz", ...
    "MaxRelativeSpacing", 0.03);

verifyEqual(testCase, result.kind, "educational_touchstone_design_frequency_grid");
verifyEqual(testCase, result.policy, "readable_touchstone_design_frequency_bracket_gate");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.indexBase, 1);
verifyEqual(testCase, result.lowerIndex, 2);
verifyEqual(testCase, result.upperIndex, 3);
verifyEqual(testCase, result.bracketGapRel, 0.02, "AbsTol", 1e-14);
verifyTrue(testCase, result.checks.designFrequencyBracketed);
verifyTrue(testCase, result.checks.designSpacingOk);

exact = educationalTouchstoneDesignFrequencyGrid([0.9; 1.0; 1.1], 1.0);
verifyEqual(testCase, exact.status, "ok");
verifyEqual(testCase, exact.lowerIndex, 2);
verifyEqual(testCase, exact.upperIndex, 2);
verifyEqual(testCase, exact.bracketGapHz, 0.0, "AbsTol", 1e-14);

coarse = educationalTouchstoneDesignFrequencyGrid([0.8 1.2], 1.0, ...
    "FrequencyUnit", "GHz", ...
    "MaxRelativeSpacing", 0.05);
verifyEqual(testCase, coarse.status, "needs_attention");
verifyTrue(testCase, coarse.checks.designFrequencyBracketed);
verifyFalse(testCase, coarse.checks.designSpacingOk);
verifyTrue(testCase, any(contains(coarse.issues, "too sparse")));

outside = educationalTouchstoneDesignFrequencyGrid([0.8 0.9], 1.0, ...
    "FrequencyUnit", "GHz");
verifyEqual(testCase, outside.status, "needs_attention");
verifyFalse(testCase, outside.checks.designFrequencyBracketed);

verifyError(testCase, @() educationalTouchstoneDesignFrequencyGrid([0.9 0.9 1.0], 1.0), ...
    "educationalTouchstoneDesignFrequencyGrid:monotonic");
end


function testTouchstoneSolverReadyPreflightBundlesRowContracts(testCase)
s21 = 0.80 * exp(-1i * deg2rad(10));

result = educationalTouchstoneSolverReadyPreflight(0.05, s21, ...
    "S12", s21, ...
    "S22", 0.05, ...
    "Frequency", 1.0, ...
    "FrequencyUnit", "GHz", ...
    "DataFormat", "MA", ...
    "Z0", 50, ...
    "ReturnLossMinDb", 20.0, ...
    "VswrMax", 1.2);

verifyEqual(testCase, result.kind, "educational_touchstone_solver_ready_preflight");
verifyEqual(testCase, result.policy, "readable_touchstone_solver_ready_row_preflight");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.frequencyHz, 1.0e9, "AbsTol", 1e-6);
verifyTrue(testCase, result.checks.frequencyRecorded);
verifyTrue(testCase, result.checks.formatRecorded);
verifyTrue(testCase, result.checks.referenceImpedanceContractOk);
verifyTrue(testCase, result.checks.sparameterPassivityOk);
verifyTrue(testCase, result.checks.sparameterReciprocityOk);
verifyTrue(testCase, result.checks.portMatchContractOk);
verifyEqual(testCase, result.portMatch.returnLossDb, 26.020599913279625, "AbsTol", 1e-12);

active = educationalTouchstoneSolverReadyPreflight(0.05, 1.2, ...
    "S12", 1.2, ...
    "S22", 0.05, ...
    "Frequency", 1.0);
verifyEqual(testCase, active.status, "needs_attention");
verifyFalse(testCase, active.checks.sparameterPassivityOk);

missingFrequency = educationalTouchstoneSolverReadyPreflight(0.05, s21, ...
    "S12", s21, ...
    "S22", 0.05);
verifyEqual(testCase, missingFrequency.status, "needs_attention");
verifyFalse(testCase, missingFrequency.checks.frequencyRecorded);
end


function path = writeFixture(testCase, text)
path = string(fullfile(tempdir, "educationalOpt_" + char(java.util.UUID.randomUUID()) + ".vol"));
fid = fopen(path, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", text);
clear cleanup
testCase.addTeardown(@() delete(path));
end


function text = tetVolText()
text = join([
    "mesh3d"
    "dimension"
    "3"
    "geomtype"
    "0"
    "facedescriptors"
    "1"
    "1 1 0 1 1"
    "surfaceelements"
    "4"
    "1 1 1 0 3 1 2 3"
    "1 1 1 0 3 1 4 2"
    "1 1 1 0 3 2 4 3"
    "1 1 1 0 3 3 4 1"
    "volumeelements"
    "1"
    "1 4 1 2 3 4"
    "points"
    "4"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "pointelements"
    "0"
    "materials"
    "1"
    "1 air"
    "bcnames"
    "1"
    "1 outer"
    "endmesh"
    ], newline);
end


function text = fourTetWithInteriorNodeVolText()
text = join([
    "mesh3d"
    "dimension"
    "3"
    "geomtype"
    "0"
    "facedescriptors"
    "1"
    "1 1 0 1 1"
    "surfaceelements"
    "4"
    "1 1 1 0 3 1 2 3"
    "1 1 1 0 3 1 4 2"
    "1 1 1 0 3 1 3 4"
    "1 1 1 0 3 2 4 3"
    "volumeelements"
    "4"
    "1 4 1 2 3 5"
    "1 4 1 4 2 5"
    "1 4 1 3 4 5"
    "1 4 2 4 3 5"
    "points"
    "5"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "0.25 0.25 0.25"
    "pointelements"
    "0"
    "materials"
    "1"
    "1 air"
    "bcnames"
    "1"
    "1 outer"
    "endmesh"
    ], newline);
end
