function tests = testFirstOrderFemBemSpaces
%TESTFIRSTORDERFEMBEMSPACES Tests for Gypsilab-style .vol first-order spaces.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testSimpleApiBuildsH1HcurlRwg(testCase)
path = writeFixture(testCase, tetVolText());

m = FemBemModel(path);

verifyEqual(testCase, m.status, "vol_ready_first_order_h1_hcurl_rwg");
verifyEqual(testCase, m.h1.basis, "P1");
verifyEqual(testCase, m.hcurl.basis, "Nedelec0");
verifyEqual(testCase, m.rwg.basis, "RWG");
verifyEqual(testCase, size(m.hcurl.edges, 1), 6);
verifyEqual(testCase, numel(m.rwg.dofEdgeIds), 6);
verifyEqual(testCase, m.rwgToHcurlEdgeIds(:), (1:6).');
end


function testH1StiffnessOnUnitTetra(testCase)
path = writeFixture(testCase, tetVolText());
m = FemBemModel(path);

[K, fem] = m.h1.stiffness();
K = full(K);
expected = (1 / 6) * [ ...
     3 -1 -1 -1
    -1  1  0  0
    -1  0  1  0
    -1  0  0  1];

verifyEqual(testCase, K, expected, "AbsTol", 1e-14);
verifyEqual(testCase, sum(K, 2), zeros(4, 1), "AbsTol", 1e-14);
verifyEqual(testCase, fem.volumes, 1/6, "AbsTol", 1e-14);
end


function testBoundaryCompactionKeepsTraceToInteriorVolumeNodes(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());

m = FemBemModel(path);

verifyEqual(testCase, size(m.mesh.vtx, 1), 5);
verifyEqual(testCase, size(m.surface.vtx, 1), 4);
verifyEqual(testCase, m.surface.volNodeIds, (1:4).');
verifyEqual(testCase, size(m.trace.matrix), [4 5]);

u = (10:10:50).';
verifyEqual(testCase, m.trace * u, u(1:4));
end


function testAssembledTraceScaffoldContainsHcurlRwgMap(testCase)
path = writeFixture(testCase, tetVolText());
m = FemBemModel(path);

m = m.assemble();

verifyEqual(testCase, m.status, "operators_ready_first_order_h1_hcurl_rwg_trace");
verifyEqual(testCase, m.operators.rwgToHcurlEdgeIds(:), (1:6).');
verifyEqual(testCase, m.operators.trace.surfaceMeshId, m.trace.surfaceMeshId);
verifyEqual(testCase, m.operators.trace.artifactId, m.trace.artifactId);
verifyEqual(testCase, m.operators.trace.operatorArtifactId, m.trace.operatorArtifactId);
verifyEqual(testCase, m.operators.trace.operatorPolicy, m.trace.operatorPolicy);
verifyEqual(testCase, m.operators.trace.outputArtifactId, m.trace.outputArtifactId);
verifyEqual(testCase, m.operators.trace.outputDigest, m.trace.outputDigest);
verifyEqual(testCase, m.operators.trace.outputPath, m.trace.outputPath);
verifyEqual(testCase, m.operators.trace.observableId, m.trace.observableId);
verifyEqual(testCase, m.operators.trace.observableFamily, m.trace.observableFamily);
verifyEqual(testCase, m.operators.trace.basisSchemaId, m.trace.basisSchemaId);
verifyEqual(testCase, m.operators.trace.assemblyRuleId, m.trace.assemblyRuleId);
verifyEqual(testCase, m.operators.trace.quadratureRuleId, m.trace.quadratureRuleId);
verifyEqual(testCase, [m.operators.trace.rowIdentity.trace_row_index].', (1:4).');
verifyEqual(testCase, [m.operators.trace.rowIdentity.fem_node_id].', (1:4).');
verifyEqual(testCase, [m.operators.trace.rowIdentity.bem_node_id].', (1:4).');
verifyEqual(testCase, m.operators.trace.rowIdentity, m.trace.rowIdentity);
verifyEqual(testCase, m.operators.trace.sourceFileId, m.trace.sourceFileId);
verifyEqual(testCase, m.surface.col, ones(4, 1));
verifyEqual(testCase, m.surface.names, repmat("outer", 4, 1));
verifyEqual(testCase, [m.surface.rowIdentity.surface_triangle_index].', (1:4).');
verifyEqual(testCase, m.surface.rowIdentity(2).surface_triangle_nodes, [1 4 2]);
verifyEqual(testCase, [m.surface.rowIdentity.boundary_number].', ones(4, 1));
verifyEqual(testCase, string({m.surface.rowIdentity.boundary_name}).', repmat("outer", 4, 1));
verifyTrue(testCase, startsWith(m.mesh.sourceFileId, "sha256:"));
verifyTrue(testCase, contains(m.trace.outputArtifactId, "h1_to_scalar_bem_trace_output_p1"));
verifyTrue(testCase, startsWith(m.trace.outputDigest, "sha256:"));
verifyTrue(testCase, startsWith(m.trace.outputPath, "memory://"));
verifyTrue(testCase, contains(m.trace.observableId, "boundary_trace_observable_p1"));
verifyEqual(testCase, m.trace.observableFamily, "fem_bem_boundary_trace");
verifyEqual(testCase, m.trace.basisSchemaId, "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1");
verifyEqual(testCase, m.trace.assemblyRuleId, "first_order_tet_h1_trace_tri_p1_bem_teaching_v1");
verifyEqual(testCase, m.trace.quadratureRuleId, "tri_p1_exact_mass_regular_kernel_teaching_v1");
verifyEqual(testCase, size(m.operators.fem.stiffness), [4 4]);
verifyEqual(testCase, size(m.operators.bem.surfaceMass), [4 4]);
end


function testTraceSurfaceMassReportKeepsBoundaryEnergyVisible(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());
m = FemBemModel(path);
u = [1; 2; -1; 0.5; 9];

report = femBemTraceMassReport(m, u);

verifyEqual(testCase, report.status, "ok");
verifyEqual(testCase, report.policy, "readable_first_order_tri_tet_trace_surface_mass_gate");
verifyNotEmpty(testCase, report.meshId);
verifyNotEmpty(testCase, report.surfaceMeshId);
verifyNotEmpty(testCase, report.traceArtifactId);
verifyTrue(testCase, contains(report.surfaceMeshId, "boundary_tri_p1"));
verifyTrue(testCase, contains(report.traceArtifactId, "h1_to_scalar_bem_trace_p1"));
verifyEqual(testCase, report.traceShape, [4 5]);
verifyEqual(testCase, report.interiorNodeIds, 5);
verifyEqual(testCase, report.surfaceAreaFromMass, 1.5 + sqrt(3) / 2, "AbsTol", 1e-14);
verifyLessThan(testCase, report.energyAbsError, 1e-14);
verifyEqual(testCase, report.interiorBoundaryAction, 0, "AbsTol", 1e-14);
verifyTrue(testCase, report.checks.surfaceMassSymmetric);
verifyTrue(testCase, report.checks.energyIdentityOk);
verifyTrue(testCase, report.checks.interiorNodesHaveZeroBoundaryAction);
verifyTrue(testCase, report.checks.surfaceMeshIdRecorded);
verifyTrue(testCase, report.checks.traceArtifactIdRecorded);
end


function testNormalFluxReportConsumesOrientationSigns(testCase)
path = writeFixture(testCase, tetVolText());
m = FemBemModel(path);

report = femBemNormalFluxSignReport(m, [1 2 3]);

verifyEqual(testCase, report.status, "ok");
verifyEqual(testCase, report.policy, "readable_fem_bem_normal_flux_orientation_gate");
verifyEqual(testCase, report.resultArtifactId, "matlab_fem_bem_normal_flux_sign_report_v1");
verifyEqual(testCase, report.expectedResultArtifactId, "matlab_fem_bem_normal_flux_sign_report_v1");
verifyTrue(testCase, startsWith(report.resultDigest, "sha256:"));
verifyEqual(testCase, report.expectedResultDigest, report.resultDigest);
verifyEqual(testCase, report.normalConvention, "outward_from_volume");
verifyEqual(testCase, report.expectedNormalConvention, "outward_from_volume");
verifyEqual(testCase, report.triangleOrientationSignsToOutward, -ones(4, 1));
verifyEqual(testCase, report.storedNormalFlux, [1.5; 1.0; -3.0; 0.5], "AbsTol", 1e-14);
verifyEqual(testCase, report.orientationCorrectedNormalFlux, [-1.5; -1.0; 3.0; -0.5], "AbsTol", 1e-14);
verifyEqual(testCase, report.outwardNormalFluxReference, [-1.5; -1.0; 3.0; -0.5], "AbsTol", 1e-14);
verifyLessThan(testCase, report.maxAbsError, 1e-14);
verifyLessThan(testCase, abs(report.closedSurfaceFluxSum), 1e-14);
verifyTrue(testCase, report.checks.normalConventionRecorded);
verifyTrue(testCase, report.checks.normalConventionMatchesExpected);
verifyTrue(testCase, report.checks.resultArtifactIdRecorded);
verifyTrue(testCase, report.checks.resultArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.resultDigestRecorded);
verifyTrue(testCase, report.checks.resultDigestMatchesExpected);
verifyTrue(testCase, report.checks.correctedFluxMatchesOutwardReference);
verifyTrue(testCase, report.checks.closedSurfaceFluxBalanceOk);

bad = m;
bad.surface.orientation.triangleOrientationSignsToOutward = ones(4, 1);
badReport = femBemNormalFluxSignReport(bad, [1 2 3]);
verifyEqual(testCase, badReport.status, "needs_attention");
verifyFalse(testCase, badReport.checks.correctedFluxMatchesOutwardReference);

wrongConvention = femBemNormalFluxSignReport( ...
    m, [1 2 3], ...
    "NormalConvention", "stored_from_vol", ...
    "ExpectedNormalConvention", "outward_from_volume");
verifyEqual(testCase, wrongConvention.status, "needs_attention");
verifyFalse(testCase, wrongConvention.checks.normalConventionMatchesExpected);
verifyTrue(testCase, wrongConvention.checks.correctedFluxMatchesOutwardReference);

staleDigest = femBemNormalFluxSignReport( ...
    m, [1 2 3], ...
    "ResultDigest", "sha256:old-normal-flux-report", ...
    "ExpectedResultDigest", report.resultDigest);
verifyEqual(testCase, staleDigest.status, "needs_attention");
verifyFalse(testCase, staleDigest.checks.resultDigestMatchesExpected);
verifyTrue(testCase, staleDigest.checks.correctedFluxMatchesOutwardReference);
end


function testCouplingManifestKeepsTraceKernelAndSpacesTogether(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());
m = FemBemModel(path);
m = m.assemble();
normalFluxArtifactId = "matlab_slot360_unit_tet_normal_flux_sign_report_v1";
normalFluxReport = femBemNormalFluxSignReport( ...
    m, [1 2 3], ...
    "ResultArtifactId", normalFluxArtifactId);
verifyEqual(testCase, normalFluxReport.status, "ok");
normalFluxDigest = normalFluxReport.resultDigest;
linearSolverReportArtifactId = "matlab_slot367_fem_bem_linear_solver_report_v1";
linearSolverReportDigest = "sha256:matlab_slot367_fem_bem_linear_solver_report_v1";
linearSolverName = "minimum_norm_pinv_rank_deficient";
linearSolverTolerance = 1e-10;
linearSolverResidualNorm = 2e-13;
linearSolverIterationCount = 1;
notebookSourceArtifactId = "matlab_slot381_fem_bem_teaching_notebook_v1";
notebookSourceDigest = "sha256:matlab-slot381-fem-bem-teaching-notebook-v1";
notebookSourcePath = "docs/fem_bem/first_order_fem_bem_teaching.ipynb";
parameterSetArtifactId = "matlab_slot388_fem_bem_parameter_set_v1";
parameterSetDigest = "sha256:matlab-slot388-fem-bem-parameter-set-v1";
parameterSetPath = "docs/fem_bem/first_order_fem_bem_parameter_set.json";
objectiveObservableId = "matlab_slot388_trace_lsq_residual_objective_v1";
objectiveObservableFamily = "fem_bem_trace_least_squares_objective";
couplingConventionSchemaId = "matlab_first_order_fem_bem_coupling_convention_v1";
postprocessRowConventionSchemaId = "matlab_fem_bem_trace_lsq_row_convention_v1";
traceBasisSchemaId = "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1";

report = femBemCouplingManifest( ...
    m, ...
    "ExpectedBoundaryNumbers", 1, ...
    "ExpectedBoundaryNames", "outer", ...
    "ExpectedTraceOperatorArtifactId", m.trace.operatorArtifactId, ...
    "ExpectedTraceOutputArtifactId", m.trace.outputArtifactId, ...
    "ExpectedTraceOutputDigest", m.trace.outputDigest, ...
    "ExpectedTraceObservableId", m.trace.observableId, ...
    "ExpectedTraceObservableFamily", "fem_bem_boundary_trace", ...
    "NormalFluxArtifactId", normalFluxArtifactId, ...
    "NormalFluxDigest", normalFluxDigest, ...
    "NormalFluxConvention", "outward_from_volume", ...
    "ExpectedNormalFluxArtifactId", normalFluxArtifactId, ...
    "ExpectedNormalFluxDigest", normalFluxDigest, ...
    "ExpectedNormalFluxConvention", "outward_from_volume", ...
    "LinearSolverReportArtifactId", linearSolverReportArtifactId, ...
    "LinearSolverReportDigest", linearSolverReportDigest, ...
    "LinearSolverName", linearSolverName, ...
    "LinearSolverTolerance", linearSolverTolerance, ...
    "LinearSolverResidualNorm", linearSolverResidualNorm, ...
    "LinearSolverIterationCount", linearSolverIterationCount, ...
    "ExpectedLinearSolverReportArtifactId", linearSolverReportArtifactId, ...
    "ExpectedLinearSolverReportDigest", linearSolverReportDigest, ...
    "ExpectedLinearSolverName", linearSolverName, ...
    "ExpectedLinearSolverTolerance", linearSolverTolerance, ...
    "ExpectedLinearSolverResidualNormMax", 1e-12, ...
    "ExpectedAssemblyRuleId", m.trace.assemblyRuleId, ...
    "ExpectedQuadratureRuleId", m.trace.quadratureRuleId, ...
    "CouplingConventionSchemaId", couplingConventionSchemaId, ...
    "ExpectedCouplingConventionSchemaId", couplingConventionSchemaId, ...
    "PostprocessRowConventionSchemaId", postprocessRowConventionSchemaId, ...
    "ExpectedPostprocessRowConventionSchemaId", postprocessRowConventionSchemaId, ...
    "ExpectedTraceBasisSchemaId", traceBasisSchemaId, ...
    "ResultArtifactId", "matlab_slot344_fem_bem_manifest_result_v1", ...
    "RunStartedAt", "2026-07-01T13:50:00+09:00", ...
    "MatlabVersion", "R2026a", ...
    "ExpectedResultArtifactId", "matlab_slot344_fem_bem_manifest_result_v1", ...
    "ExpectedMatlabVersion", "R2026a", ...
    "NotebookSourceArtifactId", notebookSourceArtifactId, ...
    "NotebookSourceDigest", notebookSourceDigest, ...
    "NotebookSourcePath", notebookSourcePath, ...
    "ExpectedNotebookSourceArtifactId", notebookSourceArtifactId, ...
    "ExpectedNotebookSourceDigest", notebookSourceDigest, ...
    "ExpectedNotebookSourcePath", notebookSourcePath, ...
    "ParameterSetArtifactId", parameterSetArtifactId, ...
    "ParameterSetDigest", parameterSetDigest, ...
    "ParameterSetPath", parameterSetPath, ...
    "ExpectedParameterSetArtifactId", parameterSetArtifactId, ...
    "ExpectedParameterSetDigest", parameterSetDigest, ...
    "ExpectedParameterSetPath", parameterSetPath, ...
    "ObjectiveObservableId", objectiveObservableId, ...
    "ObjectiveObservableFamily", objectiveObservableFamily, ...
    "ExpectedObjectiveObservableId", objectiveObservableId, ...
    "ExpectedObjectiveObservableFamily", objectiveObservableFamily, ...
    "TimingBreakdown", struct( ...
        "mesh_read_s", 0.001, ...
        "trace_assembly_s", 0.002, ...
        "manifest_build_s", 0.003, ...
        "json_write_s", 0.004), ...
    "RequireLinearSolverReport", true, ...
    "RequireNormalFluxArtifact", true, ...
    "RequireNotebookSourceArtifact", true, ...
    "RequireParameterSetArtifact", true, ...
    "RequireCouplingConventionSchema", true, ...
    "RequirePostprocessRowConventionSchema", true, ...
    "RequireTraceBasisSchema", true, ...
    "RequireTraceOutputArtifact", true);

verifyEqual(testCase, report.status, "ok");
verifyEqual(testCase, report.policy, "readable_fem_bem_coupling_manifest_gate");
verifyEqual(testCase, report.couplingKind, "h1_p1_to_scalar_bem_p1_trace");
verifyEqual(testCase, report.formulationId, "laplace_single_layer_teaching");
verifyEqual(testCase, report.bemKernelFamily, "laplace_single_layer");
verifyEqual(testCase, report.couplingConventionSchemaId, couplingConventionSchemaId);
verifyEqual(testCase, report.expectedCouplingConventionSchemaId, couplingConventionSchemaId);
verifyEqual(testCase, report.postprocessRowConventionSchemaId, postprocessRowConventionSchemaId);
verifyEqual(testCase, report.expectedPostprocessRowConventionSchemaId, postprocessRowConventionSchemaId);
verifyEqual(testCase, report.traceBasisSchemaId, traceBasisSchemaId);
verifyEqual(testCase, report.operatorTraceBasisSchemaId, traceBasisSchemaId);
verifyEqual(testCase, report.expectedTraceBasisSchemaId, traceBasisSchemaId);
verifyEqual(testCase, report.bemKernelManifestId, "laplace_single_layer_static_kernel_manifest_v1");
verifyEqual(testCase, report.bemKernelStrategy, "direct_laplace_1_over_4pi_r");
verifyEqual(testCase, report.kernelTimeConvention, "static_laplace_no_time_convention");
verifyEqual(testCase, report.assemblyRuleId, "first_order_tet_h1_trace_tri_p1_bem_teaching_v1");
verifyEqual(testCase, report.operatorTraceAssemblyRuleId, "first_order_tet_h1_trace_tri_p1_bem_teaching_v1");
verifyEqual(testCase, report.quadratureRuleId, "tri_p1_exact_mass_regular_kernel_teaching_v1");
verifyEqual(testCase, report.operatorTraceQuadratureRuleId, "tri_p1_exact_mass_regular_kernel_teaching_v1");
verifyEqual(testCase, report.volumeSpace, "H1_P1_tet");
verifyEqual(testCase, report.surfaceSpace, "scalar_P1_tri");
verifyEqual(testCase, report.boundaryNumbers, ones(4, 1));
verifyEqual(testCase, report.boundaryNames, repmat("outer", 4, 1));
verifyEqual(testCase, report.boundaryRowIdentity, m.surface.rowIdentity);
verifyEqual(testCase, report.expectedBoundaryNumbers, 1);
verifyEqual(testCase, report.expectedBoundaryNames, "outer");
verifyEqual(testCase, report.traceOperatorArtifactId, m.trace.operatorArtifactId);
verifyEqual(testCase, report.operatorTraceOperatorArtifactId, m.trace.operatorArtifactId);
verifyEqual(testCase, report.traceOperatorPolicy, "one_hot_boundary_node_injection_from_vol_node_ids");
verifyEqual(testCase, report.operatorTraceOperatorPolicy, "one_hot_boundary_node_injection_from_vol_node_ids");
verifyEqual(testCase, report.traceOutputArtifactId, m.trace.outputArtifactId);
verifyEqual(testCase, report.operatorTraceOutputArtifactId, m.trace.outputArtifactId);
verifyEqual(testCase, report.traceOutputDigest, m.trace.outputDigest);
verifyEqual(testCase, report.operatorTraceOutputDigest, m.trace.outputDigest);
verifyEqual(testCase, report.traceOutputPath, m.trace.outputPath);
verifyEqual(testCase, report.operatorTraceOutputPath, m.trace.outputPath);
verifyEqual(testCase, report.traceObservableId, m.trace.observableId);
verifyEqual(testCase, report.operatorTraceObservableId, m.trace.observableId);
verifyEqual(testCase, report.traceObservableFamily, "fem_bem_boundary_trace");
verifyEqual(testCase, report.operatorTraceObservableFamily, "fem_bem_boundary_trace");
verifyEqual(testCase, report.normalFluxArtifactId, normalFluxArtifactId);
verifyEqual(testCase, report.normalFluxDigest, normalFluxDigest);
verifyEqual(testCase, report.normalFluxConvention, "outward_from_volume");
verifyEqual(testCase, report.expectedNormalFluxArtifactId, normalFluxArtifactId);
verifyEqual(testCase, report.expectedNormalFluxDigest, normalFluxDigest);
verifyEqual(testCase, report.expectedNormalFluxConvention, "outward_from_volume");
verifyEqual(testCase, report.linearSolverReportArtifactId, linearSolverReportArtifactId);
verifyEqual(testCase, report.linearSolverReportDigest, linearSolverReportDigest);
verifyEqual(testCase, report.linearSolverName, linearSolverName);
verifyEqual(testCase, report.linearSolverTolerance, linearSolverTolerance, "AbsTol", 1e-15);
verifyEqual(testCase, report.linearSolverResidualNorm, linearSolverResidualNorm, "AbsTol", 1e-18);
verifyEqual(testCase, report.linearSolverIterationCount, linearSolverIterationCount);
verifyEqual(testCase, report.expectedLinearSolverReportArtifactId, linearSolverReportArtifactId);
verifyEqual(testCase, report.expectedLinearSolverReportDigest, linearSolverReportDigest);
verifyEqual(testCase, report.expectedLinearSolverName, linearSolverName);
verifyEqual(testCase, report.expectedLinearSolverTolerance, linearSolverTolerance, "AbsTol", 1e-15);
verifyEqual(testCase, report.expectedLinearSolverResidualNormMax, 1e-12, "AbsTol", 1e-18);
verifyEqual(testCase, report.resultArtifactId, "matlab_slot344_fem_bem_manifest_result_v1");
verifyEqual(testCase, report.expectedResultArtifactId, "matlab_slot344_fem_bem_manifest_result_v1");
verifyEqual(testCase, report.runStartedAt, "2026-07-01T13:50:00+09:00");
verifyEqual(testCase, report.matlabVersion, "R2026a");
verifyEqual(testCase, report.expectedMatlabVersion, "R2026a");
verifyEqual(testCase, report.notebookSourceArtifactId, notebookSourceArtifactId);
verifyEqual(testCase, report.notebookSourceDigest, notebookSourceDigest);
verifyEqual(testCase, report.notebookSourcePath, notebookSourcePath);
verifyEqual(testCase, report.expectedNotebookSourceArtifactId, notebookSourceArtifactId);
verifyEqual(testCase, report.expectedNotebookSourceDigest, notebookSourceDigest);
verifyEqual(testCase, report.expectedNotebookSourcePath, notebookSourcePath);
verifyEqual(testCase, report.parameterSetArtifactId, parameterSetArtifactId);
verifyEqual(testCase, report.parameterSetDigest, parameterSetDigest);
verifyEqual(testCase, report.parameterSetPath, parameterSetPath);
verifyEqual(testCase, report.expectedParameterSetArtifactId, parameterSetArtifactId);
verifyEqual(testCase, report.expectedParameterSetDigest, parameterSetDigest);
verifyEqual(testCase, report.expectedParameterSetPath, parameterSetPath);
verifyEqual(testCase, report.objectiveObservableId, objectiveObservableId);
verifyEqual(testCase, report.objectiveObservableFamily, objectiveObservableFamily);
verifyEqual(testCase, report.expectedObjectiveObservableId, objectiveObservableId);
verifyEqual(testCase, report.expectedObjectiveObservableFamily, objectiveObservableFamily);
verifyEqual(testCase, report.execution.resultArtifactId, report.resultArtifactId);
verifyEqual(testCase, report.execution.linearSolverReportArtifactId, linearSolverReportArtifactId);
verifyEqual(testCase, report.execution.linearSolverReportDigest, linearSolverReportDigest);
verifyEqual(testCase, report.execution.linearSolverName, linearSolverName);
verifyEqual(testCase, report.execution.linearSolverResidualNorm, linearSolverResidualNorm, "AbsTol", 1e-18);
verifyEqual(testCase, report.execution.matlabVersion, report.matlabVersion);
verifyEqual(testCase, report.execution.notebookSourceArtifactId, notebookSourceArtifactId);
verifyEqual(testCase, report.execution.notebookSourceDigest, notebookSourceDigest);
verifyEqual(testCase, report.execution.notebookSourcePath, notebookSourcePath);
verifyEqual(testCase, report.execution.parameterSetArtifactId, parameterSetArtifactId);
verifyEqual(testCase, report.execution.parameterSetDigest, parameterSetDigest);
verifyEqual(testCase, report.execution.parameterSetPath, parameterSetPath);
verifyEqual(testCase, report.execution.objectiveObservableId, objectiveObservableId);
verifyEqual(testCase, report.execution.objectiveObservableFamily, objectiveObservableFamily);
verifyEqual(testCase, report.optimization.parameterSetArtifactId, parameterSetArtifactId);
verifyEqual(testCase, report.optimization.parameterSetDigest, parameterSetDigest);
verifyEqual(testCase, report.optimization.parameterSetPath, parameterSetPath);
verifyEqual(testCase, report.optimization.objectiveObservableId, objectiveObservableId);
verifyEqual(testCase, report.optimization.objectiveObservableFamily, objectiveObservableFamily);
verifyEqual(testCase, report.timingBreakdownNames, ...
    ["mesh_read_s"; "trace_assembly_s"; "manifest_build_s"; "json_write_s"]);
verifyEqual(testCase, report.timingBreakdownSeconds, [0.001; 0.002; 0.003; 0.004], "AbsTol", 1e-15);
verifyEqual(testCase, report.timingTotalSeconds, 0.010, "AbsTol", 1e-15);
verifyEqual(testCase, report.requireTraceOutputArtifact, true);
verifyEqual(testCase, report.requireNormalFluxArtifact, true);
verifyEqual(testCase, report.requireLinearSolverReport, true);
verifyEqual(testCase, report.requireNotebookSourceArtifact, true);
verifyEqual(testCase, report.requireParameterSetArtifact, true);
verifyEqual(testCase, report.requirePostprocessRowConventionSchema, true);
verifyEqual(testCase, report.requireTraceBasisSchema, true);
verifyEqual(testCase, report.traceShape, [4 5]);
verifyEqual(testCase, report.boundaryNodeIds, (1:4).');
verifyEqual(testCase, [report.traceRowIdentity.trace_row_index].', (1:4).');
verifyEqual(testCase, [report.traceRowIdentity.fem_node_id].', (1:4).');
verifyEqual(testCase, [report.traceRowIdentity.bem_node_id].', (1:4).');
verifyEqual(testCase, [report.trace.trace_row_identity.fem_node_id].', (1:4).');
verifyEqual(testCase, [report.operatorTraceRowIdentity.trace_row_index].', (1:4).');
verifyEqual(testCase, [report.operatorTraceRowIdentity.fem_node_id].', (1:4).');
verifyEqual(testCase, report.operatorTraceRowIdentity, report.traceRowIdentity);
verifyEqual(testCase, [report.boundaryRowIdentity.surface_triangle_index].', (1:4).');
verifyEqual(testCase, report.boundaryRowIdentity(2).surface_triangle_nodes, [1 4 2]);
verifyEqual(testCase, [report.boundaryRowIdentity.boundary_number].', ones(4, 1));
verifyEqual(testCase, string({report.boundaryRowIdentity.boundary_name}).', repmat("outer", 4, 1));
verifyTrue(testCase, startsWith(report.sourceFileId, "sha256:"));
verifyEqual(testCase, report.trace.source_file_id, report.sourceFileId);
verifyEqual(testCase, report.trace.normal_flux_artifact_id, normalFluxArtifactId);
verifyEqual(testCase, report.trace.normal_flux_digest, normalFluxDigest);
verifyEqual(testCase, report.trace.normal_flux_convention, "outward_from_volume");
verifyEqual(testCase, report.trace.linear_solver_report_artifact_id, linearSolverReportArtifactId);
verifyEqual(testCase, report.trace.linear_solver_report_digest, linearSolverReportDigest);
verifyEqual(testCase, report.trace.linear_solver_name, linearSolverName);
verifyEqual(testCase, report.trace.linear_solver_residual_norm, linearSolverResidualNorm, "AbsTol", 1e-18);
verifyEqual(testCase, report.trace.fem_bem_postprocess_row_convention_schema_id, postprocessRowConventionSchemaId);
verifyEqual(testCase, report.trace.trace_basis_schema_id, traceBasisSchemaId);
verifyEqual(testCase, report.operatorTraceSourceFileId, report.sourceFileId);
verifyTrue(testCase, report.checks.couplingKindMatchesExpected);
verifyTrue(testCase, report.checks.formulationIdMatchesExpected);
verifyTrue(testCase, report.checks.bemKernelFamilyMatchesExpected);
verifyTrue(testCase, report.checks.bemKernelManifestIdRecorded);
verifyTrue(testCase, report.checks.bemKernelManifestIdMatchesExpected);
verifyTrue(testCase, report.checks.couplingConventionSchemaIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.couplingConventionSchemaIdMatchesExpected);
verifyTrue(testCase, report.checks.postprocessRowConventionSchemaIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.postprocessRowConventionSchemaIdMatchesExpected);
verifyTrue(testCase, report.checks.traceBasisSchemaIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.traceBasisSchemaIdMatchesExpected);
verifyTrue(testCase, report.checks.operatorTraceBasisSchemaIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceBasisSchemaIdMatchesTrace);
verifyTrue(testCase, report.checks.bemKernelStrategyRecorded);
verifyTrue(testCase, report.checks.bemKernelStrategyMatchesExpected);
verifyTrue(testCase, report.checks.kernelTimeConventionRecorded);
verifyTrue(testCase, report.checks.kernelTimeConventionMatchesExpected);
verifyTrue(testCase, report.checks.assemblyRuleIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceAssemblyRuleIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceAssemblyRuleIdMatchesTrace);
verifyTrue(testCase, report.checks.assemblyRuleIdMatchesExpected);
verifyTrue(testCase, report.checks.quadratureRuleIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceQuadratureRuleIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceQuadratureRuleIdMatchesTrace);
verifyTrue(testCase, report.checks.quadratureRuleIdMatchesExpected);
verifyTrue(testCase, report.checks.volumeSpaceMatchesExpected);
verifyTrue(testCase, report.checks.surfaceSpaceMatchesExpected);
verifyTrue(testCase, report.checks.sourceFileIdRecorded);
verifyTrue(testCase, report.checks.traceSourceFileIdRecorded);
verifyTrue(testCase, report.checks.traceSourceFileIdMatchesIdentity);
verifyTrue(testCase, report.checks.operatorTraceSourceFileIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceSourceFileIdMatchesIdentity);
verifyTrue(testCase, report.checks.traceOperatorArtifactIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceOperatorArtifactIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceOperatorArtifactIdMatchesTrace);
verifyTrue(testCase, report.checks.traceOperatorPolicyRecorded);
verifyTrue(testCase, report.checks.operatorTraceOperatorPolicyRecorded);
verifyTrue(testCase, report.checks.operatorTraceOperatorPolicyMatchesTrace);
verifyTrue(testCase, report.checks.traceOperatorPolicyMatchesExpected);
verifyTrue(testCase, report.checks.traceOperatorArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.traceOutputArtifactIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.traceOutputDigestRecordedWhenRequired);
verifyTrue(testCase, report.checks.traceOutputPathRecordedWhenRequired);
verifyTrue(testCase, report.checks.operatorTraceOutputArtifactIdMatchesTrace);
verifyTrue(testCase, report.checks.operatorTraceOutputDigestMatchesTrace);
verifyTrue(testCase, report.checks.operatorTraceOutputPathMatchesTrace);
verifyTrue(testCase, report.checks.traceOutputArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.traceOutputDigestMatchesExpected);
verifyTrue(testCase, report.checks.traceObservableIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceObservableIdRecorded);
verifyTrue(testCase, report.checks.operatorTraceObservableIdMatchesTrace);
verifyTrue(testCase, report.checks.traceObservableIdMatchesExpected);
verifyTrue(testCase, report.checks.traceObservableFamilyRecorded);
verifyTrue(testCase, report.checks.operatorTraceObservableFamilyRecorded);
verifyTrue(testCase, report.checks.operatorTraceObservableFamilyMatchesTrace);
verifyTrue(testCase, report.checks.traceObservableFamilyMatchesExpected);
verifyTrue(testCase, report.checks.normalFluxConventionRecorded);
verifyTrue(testCase, report.checks.normalFluxConventionMatchesExpected);
verifyTrue(testCase, report.checks.normalFluxArtifactIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.normalFluxDigestRecordedWhenRequired);
verifyTrue(testCase, report.checks.normalFluxArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.normalFluxDigestMatchesExpected);
verifyTrue(testCase, report.checks.linearSolverReportArtifactIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.linearSolverReportDigestRecordedWhenRequired);
verifyTrue(testCase, report.checks.linearSolverNameRecordedWhenRequired);
verifyTrue(testCase, report.checks.linearSolverToleranceRecordedWhenRequired);
verifyTrue(testCase, report.checks.linearSolverResidualNormRecordedWhenRequired);
verifyTrue(testCase, report.checks.linearSolverToleranceFinitePositiveWhenPresent);
verifyTrue(testCase, report.checks.linearSolverResidualNormFiniteNonnegativeWhenPresent);
verifyTrue(testCase, report.checks.linearSolverIterationCountNonnegativeWhenPresent);
verifyTrue(testCase, report.checks.linearSolverReportArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.linearSolverReportDigestMatchesExpected);
verifyTrue(testCase, report.checks.linearSolverNameMatchesExpected);
verifyTrue(testCase, report.checks.linearSolverToleranceMatchesExpected);
verifyTrue(testCase, report.checks.linearSolverResidualNormBelowExpectedMax);
verifyTrue(testCase, report.checks.resultArtifactIdRecorded);
verifyTrue(testCase, report.checks.resultArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.runStartedAtRecorded);
verifyTrue(testCase, report.checks.runStartedAtIsoLike);
verifyTrue(testCase, report.checks.matlabVersionRecorded);
verifyTrue(testCase, report.checks.matlabVersionMatchesExpected);
verifyTrue(testCase, report.checks.notebookSourceArtifactIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.notebookSourceDigestRecordedWhenRequired);
verifyTrue(testCase, report.checks.notebookSourcePathRecordedWhenRequired);
verifyTrue(testCase, report.checks.notebookSourceArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.notebookSourceDigestMatchesExpected);
verifyTrue(testCase, report.checks.notebookSourcePathMatchesExpected);
verifyTrue(testCase, report.checks.parameterSetArtifactIdRecordedWhenRequired);
verifyTrue(testCase, report.checks.parameterSetDigestRecordedWhenRequired);
verifyTrue(testCase, report.checks.parameterSetPathRecordedWhenRequired);
verifyTrue(testCase, report.checks.parameterSetArtifactIdMatchesExpected);
verifyTrue(testCase, report.checks.parameterSetDigestMatchesExpected);
verifyTrue(testCase, report.checks.parameterSetPathMatchesExpected);
verifyTrue(testCase, report.checks.objectiveObservableIdMatchesExpected);
verifyTrue(testCase, report.checks.objectiveObservableFamilyMatchesExpected);
verifyTrue(testCase, report.checks.timingBreakdownRecorded);
verifyTrue(testCase, report.checks.timingBreakdownHasFourItems);
verifyTrue(testCase, report.checks.timingBreakdownHasAtMostFourItems);
verifyTrue(testCase, report.checks.timingBreakdownFiniteNonnegative);
verifyTrue(testCase, report.checks.boundaryNumbersRecorded);
verifyTrue(testCase, report.checks.boundaryNamesRecorded);
verifyTrue(testCase, report.checks.boundaryRowIdentityRecorded);
verifyTrue(testCase, report.checks.boundaryRowIdentityRowIndicesMatch);
verifyTrue(testCase, report.checks.boundaryRowIdentityTrianglesMatch);
verifyTrue(testCase, report.checks.boundaryRowIdentityNumbersMatch);
verifyTrue(testCase, report.checks.boundaryRowIdentityNamesMatch);
verifyTrue(testCase, report.checks.boundaryNumbersMatchExpected);
verifyTrue(testCase, report.checks.boundaryNamesMatchExpected);
verifyTrue(testCase, report.checks.traceRowsAreOneHot);
verifyTrue(testCase, report.checks.traceRowIdentityRecorded);
verifyTrue(testCase, report.checks.traceRowIdentityRowIndicesMatch);
verifyTrue(testCase, report.checks.traceRowIdentityFemNodesMatch);
verifyTrue(testCase, report.checks.traceRowIdentityBemNodesMatch);
verifyTrue(testCase, report.checks.traceRowIdentityUnique);
verifyTrue(testCase, report.checks.traceRowIdentityMatchesTraceMatrix);
verifyTrue(testCase, report.checks.operatorTraceRowIdentityRecorded);
verifyTrue(testCase, report.checks.operatorTraceRowIdentityMatchesTrace);
verifyTrue(testCase, report.checks.operatorTraceRowIdentityRowIndicesMatch);
verifyTrue(testCase, report.checks.operatorTraceRowIdentityFemNodesMatch);
verifyTrue(testCase, report.checks.operatorTraceRowIdentityBemNodesMatch);
verifyTrue(testCase, report.checks.operatorTraceRowIdentityMatchesTraceMatrix);

wrongKernel = femBemCouplingManifest( ...
    m, ...
    "ExpectedBemKernelFamily", "helmholtz_single_layer");
verifyEqual(testCase, wrongKernel.status, "needs_attention");
verifyFalse(testCase, wrongKernel.checks.bemKernelFamilyMatchesExpected);
verifyTrue(testCase, wrongKernel.checks.traceRowsAreOneHot);

% assemblyRuleId is a TraceOperator constant now, so assembly-rule drift is
% simulated by expecting a different rule id instead of mutating the model.
staleAssemblyReport = femBemCouplingManifest( ...
    m, ...
    "ExpectedAssemblyRuleId", "remote_field_assembly_v0");
verifyEqual(testCase, staleAssemblyReport.status, "needs_attention");
verifyFalse(testCase, staleAssemblyReport.checks.assemblyRuleIdMatchesExpected);
verifyTrue(testCase, staleAssemblyReport.checks.operatorTraceAssemblyRuleIdMatchesTrace);

% quadratureRuleId is a TraceOperator constant now, so the old blanked-rule
% mutation becomes an expected-mismatch drift probe.
staleQuadratureReport = femBemCouplingManifest( ...
    m, ...
    "ExpectedQuadratureRuleId", "remote_field_quadrature_v0");
verifyEqual(testCase, staleQuadratureReport.status, "needs_attention");
verifyFalse(testCase, staleQuadratureReport.checks.quadratureRuleIdMatchesExpected);

staleResultArtifact = femBemCouplingManifest( ...
    m, ...
    "ResultArtifactId", "matlab_slot344_old_result_v0", ...
    "ExpectedResultArtifactId", "matlab_slot344_fem_bem_manifest_result_v1");
verifyEqual(testCase, staleResultArtifact.status, "needs_attention");
verifyFalse(testCase, staleResultArtifact.checks.resultArtifactIdMatchesExpected);

missingTimingItems = femBemCouplingManifest( ...
    m, ...
    "TimingBreakdown", struct("mesh_read_s", 0.001));
verifyEqual(testCase, missingTimingItems.status, "needs_attention");
verifyTrue(testCase, missingTimingItems.checks.timingBreakdownRecorded);
verifyFalse(testCase, missingTimingItems.checks.timingBreakdownHasFourItems);
verifyTrue(testCase, missingTimingItems.checks.timingBreakdownHasAtMostFourItems);

tooManyTimingItems = femBemCouplingManifest( ...
    m, ...
    "TimingBreakdown", struct( ...
        "mesh_read_s", 0.001, ...
        "trace_assembly_s", 0.002, ...
        "linear_solve_s", 0.003, ...
        "manifest_build_s", 0.004, ...
        "json_write_s", 0.005));
verifyEqual(testCase, tooManyTimingItems.status, "needs_attention");
verifyTrue(testCase, tooManyTimingItems.checks.timingBreakdownHasFourItems);
verifyFalse(testCase, tooManyTimingItems.checks.timingBreakdownHasAtMostFourItems);

badRunTimestamp = femBemCouplingManifest( ...
    m, ...
    "RunStartedAt", "not-a-date");
verifyEqual(testCase, badRunTimestamp.status, "needs_attention");
verifyTrue(testCase, badRunTimestamp.checks.runStartedAtRecorded);
verifyFalse(testCase, badRunTimestamp.checks.runStartedAtIsoLike);

staleNotebookSourceDigest = femBemCouplingManifest( ...
    m, ...
    "NotebookSourceArtifactId", notebookSourceArtifactId, ...
    "NotebookSourceDigest", "sha256:old-fem-bem-teaching-notebook", ...
    "NotebookSourcePath", notebookSourcePath, ...
    "ExpectedNotebookSourceArtifactId", notebookSourceArtifactId, ...
    "ExpectedNotebookSourceDigest", notebookSourceDigest, ...
    "ExpectedNotebookSourcePath", notebookSourcePath, ...
    "RequireNotebookSourceArtifact", true);
verifyEqual(testCase, staleNotebookSourceDigest.status, "needs_attention");
verifyTrue(testCase, staleNotebookSourceDigest.checks.notebookSourceArtifactIdMatchesExpected);
verifyFalse(testCase, staleNotebookSourceDigest.checks.notebookSourceDigestMatchesExpected);
verifyTrue(testCase, staleNotebookSourceDigest.checks.notebookSourcePathRecordedWhenRequired);

missingNotebookSourcePath = femBemCouplingManifest( ...
    m, ...
    "NotebookSourceArtifactId", notebookSourceArtifactId, ...
    "NotebookSourceDigest", notebookSourceDigest, ...
    "ExpectedNotebookSourcePath", notebookSourcePath, ...
    "RequireNotebookSourceArtifact", true);
verifyEqual(testCase, missingNotebookSourcePath.status, "needs_attention");
verifyTrue(testCase, missingNotebookSourcePath.checks.notebookSourceArtifactIdRecordedWhenRequired);
verifyTrue(testCase, missingNotebookSourcePath.checks.notebookSourceDigestRecordedWhenRequired);
verifyFalse(testCase, missingNotebookSourcePath.checks.notebookSourcePathRecordedWhenRequired);
verifyFalse(testCase, missingNotebookSourcePath.checks.notebookSourcePathMatchesExpected);

staleParameterSetDigest = femBemCouplingManifest( ...
    m, ...
    "ParameterSetArtifactId", parameterSetArtifactId, ...
    "ParameterSetDigest", "sha256:old-fem-bem-parameter-set", ...
    "ParameterSetPath", parameterSetPath, ...
    "ExpectedParameterSetArtifactId", parameterSetArtifactId, ...
    "ExpectedParameterSetDigest", parameterSetDigest, ...
    "ExpectedParameterSetPath", parameterSetPath, ...
    "RequireParameterSetArtifact", true);
verifyEqual(testCase, staleParameterSetDigest.status, "needs_attention");
verifyTrue(testCase, staleParameterSetDigest.checks.parameterSetArtifactIdMatchesExpected);
verifyFalse(testCase, staleParameterSetDigest.checks.parameterSetDigestMatchesExpected);
verifyTrue(testCase, staleParameterSetDigest.checks.parameterSetPathRecordedWhenRequired);

missingParameterSetPath = femBemCouplingManifest( ...
    m, ...
    "ParameterSetArtifactId", parameterSetArtifactId, ...
    "ParameterSetDigest", parameterSetDigest, ...
    "ExpectedParameterSetPath", parameterSetPath, ...
    "RequireParameterSetArtifact", true);
verifyEqual(testCase, missingParameterSetPath.status, "needs_attention");
verifyTrue(testCase, missingParameterSetPath.checks.parameterSetArtifactIdRecordedWhenRequired);
verifyTrue(testCase, missingParameterSetPath.checks.parameterSetDigestRecordedWhenRequired);
verifyFalse(testCase, missingParameterSetPath.checks.parameterSetPathRecordedWhenRequired);
verifyFalse(testCase, missingParameterSetPath.checks.parameterSetPathMatchesExpected);

wrongObjectiveFamily = femBemCouplingManifest( ...
    m, ...
    "ObjectiveObservableId", objectiveObservableId, ...
    "ObjectiveObservableFamily", "remote_field_map", ...
    "ExpectedObjectiveObservableId", objectiveObservableId, ...
    "ExpectedObjectiveObservableFamily", objectiveObservableFamily);
verifyEqual(testCase, wrongObjectiveFamily.status, "needs_attention");
verifyTrue(testCase, wrongObjectiveFamily.checks.objectiveObservableIdMatchesExpected);
verifyFalse(testCase, wrongObjectiveFamily.checks.objectiveObservableFamilyMatchesExpected);

staleCouplingConventionSchema = femBemCouplingManifest( ...
    m, ...
    "CouplingConventionSchemaId", "matlab_fem_bem_value_only_convention_v0", ...
    "ExpectedCouplingConventionSchemaId", couplingConventionSchemaId, ...
    "RequireCouplingConventionSchema", true);
verifyEqual(testCase, staleCouplingConventionSchema.status, "needs_attention");
verifyTrue(testCase, staleCouplingConventionSchema.checks.couplingConventionSchemaIdRecordedWhenRequired);
verifyFalse(testCase, staleCouplingConventionSchema.checks.couplingConventionSchemaIdMatchesExpected);
verifyTrue(testCase, staleCouplingConventionSchema.checks.couplingKindMatchesExpected);

missingCouplingConventionSchema = femBemCouplingManifest( ...
    m, ...
    "CouplingConventionSchemaId", "", ...
    "ExpectedCouplingConventionSchemaId", couplingConventionSchemaId, ...
    "RequireCouplingConventionSchema", true);
verifyEqual(testCase, missingCouplingConventionSchema.status, "needs_attention");
verifyFalse(testCase, missingCouplingConventionSchema.checks.couplingConventionSchemaIdRecordedWhenRequired);
verifyFalse(testCase, missingCouplingConventionSchema.checks.couplingConventionSchemaIdMatchesExpected);

stalePostprocessRowConventionSchema = femBemCouplingManifest( ...
    m, ...
    "PostprocessRowConventionSchemaId", "matlab_fem_bem_scalar_residual_row_v0", ...
    "ExpectedPostprocessRowConventionSchemaId", postprocessRowConventionSchemaId, ...
    "RequirePostprocessRowConventionSchema", true);
verifyEqual(testCase, stalePostprocessRowConventionSchema.status, "needs_attention");
verifyTrue(testCase, stalePostprocessRowConventionSchema.checks.couplingConventionSchemaIdMatchesExpected);
verifyTrue(testCase, stalePostprocessRowConventionSchema.checks.postprocessRowConventionSchemaIdRecordedWhenRequired);
verifyFalse(testCase, stalePostprocessRowConventionSchema.checks.postprocessRowConventionSchemaIdMatchesExpected);

missingPostprocessRowConventionSchema = femBemCouplingManifest( ...
    m, ...
    "PostprocessRowConventionSchemaId", "", ...
    "ExpectedPostprocessRowConventionSchemaId", postprocessRowConventionSchemaId, ...
    "RequirePostprocessRowConventionSchema", true);
verifyEqual(testCase, missingPostprocessRowConventionSchema.status, "needs_attention");
verifyFalse(testCase, missingPostprocessRowConventionSchema.checks.postprocessRowConventionSchemaIdRecordedWhenRequired);
verifyFalse(testCase, missingPostprocessRowConventionSchema.checks.postprocessRowConventionSchemaIdMatchesExpected);

% basisSchemaId is a TraceOperator constant now; a stale schema is simulated
% by overriding the recorded TraceBasisSchemaId option. The manifest also
% flags the override as diverging from the operator constant.
staleTraceBasisSchema = femBemCouplingManifest( ...
    m, ...
    "TraceBasisSchemaId", "matlab_h1_trace_basis_value_only_v0", ...
    "ExpectedTraceBasisSchemaId", traceBasisSchemaId, ...
    "RequireTraceBasisSchema", true);
verifyEqual(testCase, staleTraceBasisSchema.status, "needs_attention");
verifyTrue(testCase, staleTraceBasisSchema.checks.traceBasisSchemaIdRecordedWhenRequired);
verifyFalse(testCase, staleTraceBasisSchema.checks.traceBasisSchemaIdMatchesExpected);
verifyFalse(testCase, staleTraceBasisSchema.checks.operatorTraceBasisSchemaIdMatchesTrace);
verifyTrue(testCase, staleTraceBasisSchema.checks.traceRowsAreOneHot);

staleNormalFluxDigest = femBemCouplingManifest( ...
    m, ...
    "NormalFluxArtifactId", normalFluxArtifactId, ...
    "NormalFluxDigest", "sha256:old-normal-flux-sign-report", ...
    "NormalFluxConvention", "outward_from_volume", ...
    "ExpectedNormalFluxArtifactId", normalFluxArtifactId, ...
    "ExpectedNormalFluxDigest", normalFluxDigest, ...
    "ExpectedNormalFluxConvention", "outward_from_volume", ...
    "RequireNormalFluxArtifact", true);
verifyEqual(testCase, staleNormalFluxDigest.status, "needs_attention");
verifyTrue(testCase, staleNormalFluxDigest.checks.normalFluxArtifactIdMatchesExpected);
verifyFalse(testCase, staleNormalFluxDigest.checks.normalFluxDigestMatchesExpected);
verifyTrue(testCase, staleNormalFluxDigest.checks.traceRowsAreOneHot);

missingNormalFluxArtifact = femBemCouplingManifest( ...
    m, ...
    "RequireNormalFluxArtifact", true);
verifyEqual(testCase, missingNormalFluxArtifact.status, "needs_attention");
verifyFalse(testCase, missingNormalFluxArtifact.checks.normalFluxArtifactIdRecordedWhenRequired);
verifyFalse(testCase, missingNormalFluxArtifact.checks.normalFluxDigestRecordedWhenRequired);
verifyTrue(testCase, missingNormalFluxArtifact.checks.traceRowsAreOneHot);

wrongNormalFluxConvention = femBemCouplingManifest( ...
    m, ...
    "NormalFluxArtifactId", normalFluxArtifactId, ...
    "NormalFluxDigest", normalFluxDigest, ...
    "NormalFluxConvention", "stored_triangle_orientation", ...
    "ExpectedNormalFluxConvention", "outward_from_volume", ...
    "RequireNormalFluxArtifact", true);
verifyEqual(testCase, wrongNormalFluxConvention.status, "needs_attention");
verifyFalse(testCase, wrongNormalFluxConvention.checks.normalFluxConventionMatchesExpected);
verifyTrue(testCase, wrongNormalFluxConvention.checks.normalFluxDigestRecordedWhenRequired);

staleSolverReportDigest = femBemCouplingManifest( ...
    m, ...
    "LinearSolverReportArtifactId", linearSolverReportArtifactId, ...
    "LinearSolverReportDigest", "sha256:old-linear-solver-report", ...
    "LinearSolverName", linearSolverName, ...
    "LinearSolverTolerance", linearSolverTolerance, ...
    "LinearSolverResidualNorm", linearSolverResidualNorm, ...
    "LinearSolverIterationCount", linearSolverIterationCount, ...
    "ExpectedLinearSolverReportArtifactId", linearSolverReportArtifactId, ...
    "ExpectedLinearSolverReportDigest", linearSolverReportDigest, ...
    "ExpectedLinearSolverName", linearSolverName, ...
    "ExpectedLinearSolverTolerance", linearSolverTolerance, ...
    "ExpectedLinearSolverResidualNormMax", 1e-12, ...
    "RequireLinearSolverReport", true);
verifyEqual(testCase, staleSolverReportDigest.status, "needs_attention");
verifyTrue(testCase, staleSolverReportDigest.checks.linearSolverReportArtifactIdMatchesExpected);
verifyFalse(testCase, staleSolverReportDigest.checks.linearSolverReportDigestMatchesExpected);
verifyTrue(testCase, staleSolverReportDigest.checks.linearSolverResidualNormBelowExpectedMax);

missingSolverReport = femBemCouplingManifest( ...
    m, ...
    "RequireLinearSolverReport", true);
verifyEqual(testCase, missingSolverReport.status, "needs_attention");
verifyFalse(testCase, missingSolverReport.checks.linearSolverReportArtifactIdRecordedWhenRequired);
verifyFalse(testCase, missingSolverReport.checks.linearSolverReportDigestRecordedWhenRequired);
verifyFalse(testCase, missingSolverReport.checks.linearSolverResidualNormRecordedWhenRequired);
verifyTrue(testCase, missingSolverReport.checks.traceRowsAreOneHot);

badSolverResidual = femBemCouplingManifest( ...
    m, ...
    "LinearSolverReportArtifactId", linearSolverReportArtifactId, ...
    "LinearSolverReportDigest", linearSolverReportDigest, ...
    "LinearSolverName", linearSolverName, ...
    "LinearSolverTolerance", linearSolverTolerance, ...
    "LinearSolverResidualNorm", 2e-6, ...
    "LinearSolverIterationCount", linearSolverIterationCount, ...
    "ExpectedLinearSolverResidualNormMax", 1e-12, ...
    "RequireLinearSolverReport", true);
verifyEqual(testCase, badSolverResidual.status, "needs_attention");
verifyFalse(testCase, badSolverResidual.checks.linearSolverResidualNormBelowExpectedMax);
verifyTrue(testCase, badSolverResidual.checks.linearSolverReportDigestRecordedWhenRequired);

wrongSolverName = femBemCouplingManifest( ...
    m, ...
    "LinearSolverReportArtifactId", linearSolverReportArtifactId, ...
    "LinearSolverReportDigest", linearSolverReportDigest, ...
    "LinearSolverName", "direct_lu", ...
    "LinearSolverTolerance", linearSolverTolerance, ...
    "LinearSolverResidualNorm", linearSolverResidualNorm, ...
    "LinearSolverIterationCount", linearSolverIterationCount, ...
    "ExpectedLinearSolverName", linearSolverName, ...
    "RequireLinearSolverReport", true);
verifyEqual(testCase, wrongSolverName.status, "needs_attention");
verifyFalse(testCase, wrongSolverName.checks.linearSolverNameMatchesExpected);
verifyTrue(testCase, wrongSolverName.checks.linearSolverResidualNormFiniteNonnegativeWhenPresent);

lowFrequencyKernel = femBemCouplingManifest( ...
    m, ...
    "FormulationId", "acoustic_low_frequency_single_layer_teaching", ...
    "BemKernelFamily", "helmholtz_single_layer", ...
    "BemKernelManifestId", "lf_helmholtz_k1e-9_expm1_taylor_v1", ...
    "BemKernelStrategy", "laplace_plus_expm1_taylor_correction", ...
    "KernelTimeConvention", "exp(+i*k*r) MATLAB teaching convention", ...
    "ExpectedFormulationId", "acoustic_low_frequency_single_layer_teaching", ...
    "ExpectedBemKernelFamily", "helmholtz_single_layer", ...
    "ExpectedBemKernelManifestId", "lf_helmholtz_k1e-9_expm1_taylor_v1", ...
    "ExpectedBemKernelStrategy", "laplace_plus_expm1_taylor_correction", ...
    "ExpectedKernelTimeConvention", "exp(+i*k*r) MATLAB teaching convention");
verifyEqual(testCase, lowFrequencyKernel.status, "ok");
verifyEqual(testCase, lowFrequencyKernel.bemKernelManifestId, ...
    "lf_helmholtz_k1e-9_expm1_taylor_v1");
verifyTrue(testCase, lowFrequencyKernel.checks.bemKernelManifestIdMatchesExpected);
verifyTrue(testCase, lowFrequencyKernel.checks.bemKernelStrategyMatchesExpected);
verifyTrue(testCase, lowFrequencyKernel.checks.kernelTimeConventionMatchesExpected);

staleKernelManifest = femBemCouplingManifest( ...
    m, ...
    "BemKernelFamily", "helmholtz_single_layer", ...
    "BemKernelManifestId", "old_direct_helmholtz_manifest", ...
    "BemKernelStrategy", "laplace_plus_expm1_taylor_correction", ...
    "KernelTimeConvention", "exp(+i*k*r) MATLAB teaching convention", ...
    "ExpectedBemKernelFamily", "helmholtz_single_layer", ...
    "ExpectedBemKernelManifestId", "lf_helmholtz_k1e-9_expm1_taylor_v1", ...
    "ExpectedBemKernelStrategy", "laplace_plus_expm1_taylor_correction", ...
    "ExpectedKernelTimeConvention", "exp(+i*k*r) MATLAB teaching convention");
verifyEqual(testCase, staleKernelManifest.status, "needs_attention");
verifyFalse(testCase, staleKernelManifest.checks.bemKernelManifestIdMatchesExpected);
verifyTrue(testCase, staleKernelManifest.checks.bemKernelStrategyMatchesExpected);

wrongKernelStrategy = femBemCouplingManifest( ...
    m, ...
    "BemKernelFamily", "helmholtz_single_layer", ...
    "BemKernelManifestId", "lf_helmholtz_k1e-9_expm1_taylor_v1", ...
    "BemKernelStrategy", "direct_exp_minus_laplace", ...
    "KernelTimeConvention", "exp(+i*k*r) MATLAB teaching convention", ...
    "ExpectedBemKernelFamily", "helmholtz_single_layer", ...
    "ExpectedBemKernelManifestId", "lf_helmholtz_k1e-9_expm1_taylor_v1", ...
    "ExpectedBemKernelStrategy", "laplace_plus_expm1_taylor_correction");
verifyEqual(testCase, wrongKernelStrategy.status, "needs_attention");
verifyTrue(testCase, wrongKernelStrategy.checks.bemKernelManifestIdMatchesExpected);
verifyFalse(testCase, wrongKernelStrategy.checks.bemKernelStrategyMatchesExpected);

badSource = m.assemble();
badSource.trace.sourceFileId = "sha256:stale_vol_source";
badSourceReport = femBemCouplingManifest(badSource);
verifyEqual(testCase, badSourceReport.status, "needs_attention");
verifyFalse(testCase, badSourceReport.checks.traceSourceFileIdMatchesIdentity);
verifyTrue(testCase, badSourceReport.checks.operatorTraceSourceFileIdMatchesIdentity);

badBoundaryName = m.assemble();
badBoundaryName.surface.names(2) = "coil";
badBoundaryNameReport = femBemCouplingManifest( ...
    badBoundaryName, ...
    "ExpectedBoundaryNames", "outer");
verifyEqual(testCase, badBoundaryNameReport.status, "needs_attention");
verifyFalse(testCase, badBoundaryNameReport.checks.boundaryNamesMatchExpected);

badBoundaryNumber = m.assemble();
badBoundaryNumber.surface.col(2) = 2;
badBoundaryNumberReport = femBemCouplingManifest( ...
    badBoundaryNumber, ...
    "ExpectedBoundaryNumbers", 1);
verifyEqual(testCase, badBoundaryNumberReport.status, "needs_attention");
verifyFalse(testCase, badBoundaryNumberReport.checks.boundaryRowIdentityNumbersMatch);
verifyFalse(testCase, badBoundaryNumberReport.checks.boundaryNumbersMatchExpected);

% Boundary row identity is single-sourced on SurfaceMesh now, so the drift
% probe mutates the surface row identity and the single-source check flags it.
badBoundaryRowIdentity = m.assemble();
badBoundaryRowIdentity.surface.rowIdentity(2).surface_triangle_nodes = [1 4 3];
badBoundaryRowIdentityReport = femBemCouplingManifest(badBoundaryRowIdentity);
verifyEqual(testCase, badBoundaryRowIdentityReport.status, "needs_attention");
verifyFalse(testCase, badBoundaryRowIdentityReport.checks.boundaryRowIdentityTrianglesMatch);

badTraceOperator = m.assemble();
badTraceOperator.operators.trace.operatorArtifactId = ...
    "netgen_vol:stale:h1_to_scalar_bem_trace_operator_p1";
badTraceOperatorReport = femBemCouplingManifest(badTraceOperator);
verifyEqual(testCase, badTraceOperatorReport.status, "needs_attention");
verifyFalse(testCase, badTraceOperatorReport.checks.operatorTraceOperatorArtifactIdMatchesTrace);

% operatorPolicy is a TraceOperator constant now; policy drift is simulated
% by expecting a different policy id instead of mutating the model.
badTracePolicyReport = femBemCouplingManifest( ...
    m, ...
    "ExpectedTraceOperatorPolicy", "remote_field_observation_map");
verifyEqual(testCase, badTracePolicyReport.status, "needs_attention");
verifyFalse(testCase, badTracePolicyReport.checks.traceOperatorPolicyMatchesExpected);

badTraceOutput = m.assemble();
badTraceOutput.operators.trace.outputArtifactId = ...
    "netgen_vol:stale:h1_to_scalar_bem_trace_output_p1";
badTraceOutputReport = femBemCouplingManifest( ...
    badTraceOutput, ...
    "ExpectedTraceOutputArtifactId", m.trace.outputArtifactId);
verifyEqual(testCase, badTraceOutputReport.status, "needs_attention");
verifyFalse(testCase, badTraceOutputReport.checks.operatorTraceOutputArtifactIdMatchesTrace);
verifyTrue(testCase, badTraceOutputReport.checks.traceOutputArtifactIdMatchesExpected);

badTraceOutputDigest = m.assemble();
badTraceOutputDigestReport = femBemCouplingManifest( ...
    badTraceOutputDigest, ...
    "ExpectedTraceOutputDigest", "sha256:stale_trace_output_digest");
verifyEqual(testCase, badTraceOutputDigestReport.status, "needs_attention");
verifyFalse(testCase, badTraceOutputDigestReport.checks.traceOutputDigestMatchesExpected);

missingTraceOutputPath = m.assemble();
missingTraceOutputPath.trace.outputPath = "";
missingTraceOutputPath.operators.trace.outputPath = "";
missingTraceOutputPathReport = femBemCouplingManifest( ...
    missingTraceOutputPath, ...
    "RequireTraceOutputArtifact", true);
verifyEqual(testCase, missingTraceOutputPathReport.status, "needs_attention");
verifyFalse(testCase, missingTraceOutputPathReport.checks.traceOutputPathRecordedWhenRequired);

badTraceObservable = m.assemble();
badTraceObservable.operators.trace.observableId = ...
    "stale_remote_field_map_observable";
badTraceObservableReport = femBemCouplingManifest( ...
    badTraceObservable, ...
    "ExpectedTraceObservableId", m.trace.observableId);
verifyEqual(testCase, badTraceObservableReport.status, "needs_attention");
verifyFalse(testCase, badTraceObservableReport.checks.operatorTraceObservableIdMatchesTrace);
verifyTrue(testCase, badTraceObservableReport.checks.traceObservableIdMatchesExpected);

% observableFamily is a TraceOperator constant now; family drift is simulated
% by expecting a different family instead of mutating the model.
badTraceObservableFamilyReport = femBemCouplingManifest( ...
    m, ...
    "ExpectedTraceObservableFamily", "remote_field_map");
verifyEqual(testCase, badTraceObservableFamilyReport.status, "needs_attention");
verifyTrue(testCase, badTraceObservableFamilyReport.checks.operatorTraceObservableFamilyMatchesTrace);
verifyFalse(testCase, badTraceObservableFamilyReport.checks.traceObservableFamilyMatchesExpected);

badTrace = m.assemble();
badTrace.trace.femNodeIds(3) = 4;
badReport = femBemCouplingManifest(badTrace);
verifyEqual(testCase, badReport.status, "needs_attention");
verifyFalse(testCase, badReport.checks.traceRowIdentityUnique);
verifyFalse(testCase, badReport.checks.traceRowIdentityMatchesTraceMatrix);

badOperatorRowIdentity = m.assemble();
badOperatorRowIdentity.operators.trace.rowIdentity(3).fem_node_id = 4;
badOperatorRowIdentityReport = femBemCouplingManifest(badOperatorRowIdentity);
verifyEqual(testCase, badOperatorRowIdentityReport.status, "needs_attention");
verifyFalse(testCase, badOperatorRowIdentityReport.checks.operatorTraceRowIdentityMatchesTrace);
verifyFalse(testCase, badOperatorRowIdentityReport.checks.operatorTraceRowIdentityFemNodesMatch);
verifyFalse(testCase, badOperatorRowIdentityReport.checks.operatorTraceRowIdentityMatchesTraceMatrix);
end


function path = writeFixture(testCase, text)
path = string(fullfile(tempdir, "firstOrderFemBem_" + char(java.util.UUID.randomUUID()) + ".vol"));
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
