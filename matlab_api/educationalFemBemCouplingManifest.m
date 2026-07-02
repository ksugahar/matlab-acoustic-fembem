function report = educationalFemBemCouplingManifest(model, options)
%EDUCATIONALFEMBEMCOUPLINGMANIFEST Package readable FEM/BEM coupling identity.
%
% This helper stops before solving.  It makes the trace handoff explicit so
% students can see which volume space, surface space, formulation, and BEM
% kernel family belong to the same .vol-derived FEM/BEM package.

arguments
    model struct
    options.CouplingKind (1,1) string = "h1_p1_to_scalar_bem_p1_trace"
    options.FormulationId (1,1) string = "laplace_single_layer_teaching"
    options.BemKernelFamily (1,1) string = "laplace_single_layer"
    options.CouplingConventionSchemaId (1,1) string = "matlab_first_order_fem_bem_coupling_convention_v1"
    options.PostprocessRowConventionSchemaId (1,1) string = "matlab_fem_bem_trace_lsq_row_convention_v1"
    options.TraceBasisSchemaId (1,1) string = ""
    options.BemKernelManifestId (1,1) string = "laplace_single_layer_static_kernel_manifest_v1"
    options.BemKernelStrategy (1,1) string = "direct_laplace_1_over_4pi_r"
    options.KernelTimeConvention (1,1) string = "static_laplace_no_time_convention"
    options.AssemblyRuleId (1,1) string = "first_order_tet_h1_trace_tri_p1_bem_teaching_v1"
    options.QuadratureRuleId (1,1) string = "tri_p1_exact_mass_regular_kernel_teaching_v1"
    options.VolumeSpace (1,1) string = "H1_P1_tet"
    options.SurfaceSpace (1,1) string = "scalar_P1_tri"
    options.ExpectedCouplingKind (1,1) string = ""
    options.ExpectedFormulationId (1,1) string = ""
    options.ExpectedBemKernelFamily (1,1) string = ""
    options.ExpectedCouplingConventionSchemaId (1,1) string = ""
    options.ExpectedPostprocessRowConventionSchemaId (1,1) string = ""
    options.ExpectedTraceBasisSchemaId (1,1) string = ""
    options.ExpectedBemKernelManifestId (1,1) string = ""
    options.ExpectedBemKernelStrategy (1,1) string = ""
    options.ExpectedKernelTimeConvention (1,1) string = ""
    options.ExpectedAssemblyRuleId (1,1) string = ""
    options.ExpectedQuadratureRuleId (1,1) string = ""
    options.ExpectedVolumeSpace (1,1) string = ""
    options.ExpectedSurfaceSpace (1,1) string = ""
    options.ExpectedBoundaryNumbers (:,1) double = []
    options.ExpectedBoundaryNames (:,1) string = strings(0, 1)
    options.ExpectedTraceOperatorArtifactId (1,1) string = ""
    options.ExpectedTraceOperatorPolicy (1,1) string = ""
    options.ExpectedTraceOutputArtifactId (1,1) string = ""
    options.ExpectedTraceOutputDigest (1,1) string = ""
    options.ExpectedTraceObservableId (1,1) string = ""
    options.ExpectedTraceObservableFamily (1,1) string = ""
    options.NormalFluxArtifactId (1,1) string = ""
    options.NormalFluxDigest (1,1) string = ""
    options.NormalFluxConvention (1,1) string = "outward_from_volume"
    options.ExpectedNormalFluxArtifactId (1,1) string = ""
    options.ExpectedNormalFluxDigest (1,1) string = ""
    options.ExpectedNormalFluxConvention (1,1) string = ""
    options.LinearSolverReportArtifactId (1,1) string = ""
    options.LinearSolverReportDigest (1,1) string = ""
    options.LinearSolverName (1,1) string = ""
    options.LinearSolverTolerance (1,1) double = NaN
    options.LinearSolverResidualNorm (1,1) double = NaN
    options.LinearSolverIterationCount (1,1) double = NaN
    options.ExpectedLinearSolverReportArtifactId (1,1) string = ""
    options.ExpectedLinearSolverReportDigest (1,1) string = ""
    options.ExpectedLinearSolverName (1,1) string = ""
    options.ExpectedLinearSolverTolerance (1,1) double = NaN
    options.ExpectedLinearSolverResidualNormMax (1,1) double = NaN
    options.ResultArtifactId (1,1) string = "matlab_fem_bem_coupling_manifest_result_v1"
    options.RunStartedAt (1,1) string = ""
    options.MatlabVersion (1,1) string = ""
    options.TimingBreakdown (1,1) struct = struct( ...
        "mesh_read_s", 0, ...
        "trace_assembly_s", 0, ...
        "manifest_build_s", 0, ...
        "json_write_s", 0)
    options.ExpectedResultArtifactId (1,1) string = ""
    options.ExpectedMatlabVersion (1,1) string = ""
    options.NotebookSourceArtifactId (1,1) string = ""
    options.NotebookSourceDigest (1,1) string = ""
    options.NotebookSourcePath (1,1) string = ""
    options.ExpectedNotebookSourceArtifactId (1,1) string = ""
    options.ExpectedNotebookSourceDigest (1,1) string = ""
    options.ExpectedNotebookSourcePath (1,1) string = ""
    options.ParameterSetArtifactId (1,1) string = ""
    options.ParameterSetDigest (1,1) string = ""
    options.ParameterSetPath (1,1) string = ""
    options.ExpectedParameterSetArtifactId (1,1) string = ""
    options.ExpectedParameterSetDigest (1,1) string = ""
    options.ExpectedParameterSetPath (1,1) string = ""
    options.ObjectiveObservableId (1,1) string = ""
    options.ObjectiveObservableFamily (1,1) string = ""
    options.ExpectedObjectiveObservableId (1,1) string = ""
    options.ExpectedObjectiveObservableFamily (1,1) string = ""
    options.RequireTraceOutputArtifact (1,1) logical = false
    options.RequireNormalFluxArtifact (1,1) logical = false
    options.RequireLinearSolverReport (1,1) logical = false
    options.RequireNotebookSourceArtifact (1,1) logical = false
    options.RequireParameterSetArtifact (1,1) logical = false
    options.RequireCouplingConventionSchema (1,1) logical = false
    options.RequirePostprocessRowConventionSchema (1,1) logical = false
    options.RequireTraceBasisSchema (1,1) logical = false
end

if ~isfield(model, "operators") || ~isfield(model.operators, "trace")
    model = assembleFirstOrderFemBemTrace(model);
end

expectedCouplingKind = defaultExpected(options.ExpectedCouplingKind, options.CouplingKind);
expectedFormulationId = defaultExpected(options.ExpectedFormulationId, options.FormulationId);
expectedBemKernelFamily = defaultExpected(options.ExpectedBemKernelFamily, options.BemKernelFamily);
expectedCouplingConventionSchemaId = defaultExpected( ...
    options.ExpectedCouplingConventionSchemaId, ...
    options.CouplingConventionSchemaId);
expectedPostprocessRowConventionSchemaId = defaultExpected( ...
    options.ExpectedPostprocessRowConventionSchemaId, ...
    options.PostprocessRowConventionSchemaId);
traceBasisSchemaId = defaultExpected( ...
    options.TraceBasisSchemaId, ...
    getStringField(model.trace, "traceBasisSchemaId"));
expectedTraceBasisSchemaId = defaultExpected( ...
    options.ExpectedTraceBasisSchemaId, ...
    traceBasisSchemaId);
expectedBemKernelManifestId = defaultExpected( ...
    options.ExpectedBemKernelManifestId, ...
    options.BemKernelManifestId);
expectedBemKernelStrategy = defaultExpected( ...
    options.ExpectedBemKernelStrategy, ...
    options.BemKernelStrategy);
expectedKernelTimeConvention = defaultExpected( ...
    options.ExpectedKernelTimeConvention, ...
    options.KernelTimeConvention);
expectedAssemblyRuleId = defaultExpected( ...
    options.ExpectedAssemblyRuleId, ...
    options.AssemblyRuleId);
expectedQuadratureRuleId = defaultExpected( ...
    options.ExpectedQuadratureRuleId, ...
    options.QuadratureRuleId);
expectedVolumeSpace = defaultExpected(options.ExpectedVolumeSpace, options.VolumeSpace);
expectedSurfaceSpace = defaultExpected(options.ExpectedSurfaceSpace, options.SurfaceSpace);
expectedTraceOperatorArtifactId = options.ExpectedTraceOperatorArtifactId;
expectedTraceOperatorPolicy = defaultExpected( ...
    options.ExpectedTraceOperatorPolicy, ...
    "one_hot_boundary_node_injection_from_vol_node_ids");
expectedTraceObservableId = options.ExpectedTraceObservableId;
expectedTraceObservableFamily = defaultExpected( ...
    options.ExpectedTraceObservableFamily, ...
    "fem_bem_boundary_trace");
expectedNormalFluxConvention = defaultExpected( ...
    options.ExpectedNormalFluxConvention, ...
    options.NormalFluxConvention);
expectedLinearSolverReportArtifactId = defaultExpected( ...
    options.ExpectedLinearSolverReportArtifactId, ...
    options.LinearSolverReportArtifactId);
expectedLinearSolverReportDigest = defaultExpected( ...
    options.ExpectedLinearSolverReportDigest, ...
    options.LinearSolverReportDigest);
expectedLinearSolverName = defaultExpected( ...
    options.ExpectedLinearSolverName, ...
    options.LinearSolverName);
expectedLinearSolverTolerance = defaultExpectedDouble( ...
    options.ExpectedLinearSolverTolerance, ...
    options.LinearSolverTolerance);
expectedResultArtifactId = defaultExpected( ...
    options.ExpectedResultArtifactId, ...
    options.ResultArtifactId);
expectedMatlabVersion = options.ExpectedMatlabVersion;
expectedNotebookSourceArtifactId = options.ExpectedNotebookSourceArtifactId;
expectedNotebookSourceDigest = options.ExpectedNotebookSourceDigest;
expectedNotebookSourcePath = options.ExpectedNotebookSourcePath;
expectedParameterSetArtifactId = options.ExpectedParameterSetArtifactId;
expectedParameterSetDigest = options.ExpectedParameterSetDigest;
expectedParameterSetPath = options.ExpectedParameterSetPath;
expectedObjectiveObservableId = options.ExpectedObjectiveObservableId;
expectedObjectiveObservableFamily = options.ExpectedObjectiveObservableFamily;
runStartedAt = options.RunStartedAt;
if strlength(runStartedAt) == 0
    runStartedAt = string(datetime("now", "TimeZone", "local", ...
        "Format", "yyyy-MM-dd'T'HH:mm:ssXXX"));
end
matlabVersion = options.MatlabVersion;
if strlength(matlabVersion) == 0
    matlabVersion = string(version);
end
[timingNames, timingValues] = timingBreakdownArrays(options.TimingBreakdown);

T = full(model.operators.trace.matrix);
traceNodeIds = model.operators.trace.femNodeIds(:);
surfaceTriangles = model.trace.surfaceTrianglesGlobal;
identitySourceFileId = getStringField(model.identity, "sourceFileId");
traceSourceFileId = getStringField(model.trace, "sourceFileId");
operatorTraceSourceFileId = getStringField(model.operators.trace, "sourceFileId");
traceOperatorArtifactId = getStringField(model.trace, "traceOperatorArtifactId");
operatorTraceOperatorArtifactId = getStringField(model.operators.trace, "traceOperatorArtifactId");
traceOperatorPolicy = getStringField(model.trace, "traceOperatorPolicy");
operatorTraceOperatorPolicy = getStringField(model.operators.trace, "traceOperatorPolicy");
traceOutputArtifactId = getStringField(model.trace, "traceOutputArtifactId");
operatorTraceOutputArtifactId = getStringField(model.operators.trace, "traceOutputArtifactId");
traceOutputDigest = getStringField(model.trace, "traceOutputDigest");
operatorTraceOutputDigest = getStringField(model.operators.trace, "traceOutputDigest");
traceOutputPath = getStringField(model.trace, "traceOutputPath");
operatorTraceOutputPath = getStringField(model.operators.trace, "traceOutputPath");
traceObservableId = getStringField(model.trace, "traceObservableId");
operatorTraceObservableId = getStringField(model.operators.trace, "traceObservableId");
traceObservableFamily = getStringField(model.trace, "traceObservableFamily");
operatorTraceObservableFamily = getStringField(model.operators.trace, "traceObservableFamily");
operatorTraceBasisSchemaId = getStringField(model.operators.trace, "traceBasisSchemaId");
traceAssemblyRuleId = getStringField(model.trace, "assemblyRuleId");
operatorTraceAssemblyRuleId = getStringField(model.operators.trace, "assemblyRuleId");
traceQuadratureRuleId = getStringField(model.trace, "quadratureRuleId");
operatorTraceQuadratureRuleId = getStringField(model.operators.trace, "quadratureRuleId");
traceRowIdentity = makeTraceRowIdentity(traceNodeIds);
operatorTraceRowIdentity = getStructVectorField(model.operators.trace, "traceRowIdentity");
rowNonzeroCount = sum(abs(T) > 0, 2);
rowMax = zeros(size(T, 1), 1);
traceMatrixNodeIds = zeros(size(T, 1), 1);
for k = 1:size(T, 1)
    rowMax(k) = max(abs(T(k, :)));
    nonzero = find(abs(T(k, :)) > 0);
    if isscalar(nonzero)
        traceMatrixNodeIds(k) = nonzero;
    else
        traceMatrixNodeIds(k) = NaN;
    end
end
traceRowIndex = [traceRowIdentity.trace_row_index].';
traceRowFemNodeIds = [traceRowIdentity.fem_node_id].';
traceRowBemNodeIds = [traceRowIdentity.bem_node_id].';
operatorTraceRowIndex = [operatorTraceRowIdentity.trace_row_index].';
operatorTraceRowFemNodeIds = [operatorTraceRowIdentity.fem_node_id].';
operatorTraceRowBemNodeIds = [operatorTraceRowIdentity.bem_node_id].';
boundaryNumbers = getNumericVectorField(model.trace, "boundaryNumbers");
boundaryNames = getStringVectorField(model.trace, "boundaryNames");
operatorBoundaryNumbers = getNumericVectorField(model.operators.trace, "boundaryNumbers");
operatorBoundaryNames = getStringVectorField(model.operators.trace, "boundaryNames");
boundaryRowIdentity = getBoundaryRowIdentityField(model.trace, "boundaryRowIdentity");
operatorBoundaryRowIdentity = getBoundaryRowIdentityField(model.operators.trace, "boundaryRowIdentity");
[boundaryRowIndices, boundaryRowTriangles, boundaryRowNumbers, boundaryRowNames] = ...
    boundaryRowIdentityArrays(boundaryRowIdentity);
[operatorBoundaryRowIndices, operatorBoundaryRowTriangles, operatorBoundaryRowNumbers, operatorBoundaryRowNames] = ...
    boundaryRowIdentityArrays(operatorBoundaryRowIdentity);
expectedBoundaryNumbers = sort(unique(options.ExpectedBoundaryNumbers(:)));
expectedBoundaryNames = sort(unique(options.ExpectedBoundaryNames(:)));

report = struct();
report.policy = "readable_fem_bem_coupling_manifest_gate";
report.meshId = model.identity.meshId;
report.exportId = model.identity.meshId;
report.sourcePath = model.identity.sourcePath;
report.sourceFileId = identitySourceFileId;
report.sourceFormat = model.identity.sourceFormat;
report.surfaceMeshId = model.trace.surfaceMeshId;
report.traceArtifactId = model.trace.traceArtifactId;
report.traceOperatorArtifactId = traceOperatorArtifactId;
report.operatorTraceOperatorArtifactId = operatorTraceOperatorArtifactId;
report.traceOperatorPolicy = traceOperatorPolicy;
report.operatorTraceOperatorPolicy = operatorTraceOperatorPolicy;
report.expectedTraceOperatorArtifactId = expectedTraceOperatorArtifactId;
report.expectedTraceOperatorPolicy = expectedTraceOperatorPolicy;
report.traceOutputArtifactId = traceOutputArtifactId;
report.operatorTraceOutputArtifactId = operatorTraceOutputArtifactId;
report.traceOutputDigest = traceOutputDigest;
report.operatorTraceOutputDigest = operatorTraceOutputDigest;
report.traceOutputPath = traceOutputPath;
report.operatorTraceOutputPath = operatorTraceOutputPath;
report.expectedTraceOutputArtifactId = options.ExpectedTraceOutputArtifactId;
report.expectedTraceOutputDigest = options.ExpectedTraceOutputDigest;
report.traceObservableId = traceObservableId;
report.operatorTraceObservableId = operatorTraceObservableId;
report.expectedTraceObservableId = expectedTraceObservableId;
report.traceObservableFamily = traceObservableFamily;
report.operatorTraceObservableFamily = operatorTraceObservableFamily;
report.expectedTraceObservableFamily = expectedTraceObservableFamily;
report.normalFluxArtifactId = options.NormalFluxArtifactId;
report.normalFluxDigest = options.NormalFluxDigest;
report.normalFluxConvention = options.NormalFluxConvention;
report.expectedNormalFluxArtifactId = options.ExpectedNormalFluxArtifactId;
report.expectedNormalFluxDigest = options.ExpectedNormalFluxDigest;
report.expectedNormalFluxConvention = expectedNormalFluxConvention;
report.linearSolverReportArtifactId = options.LinearSolverReportArtifactId;
report.linearSolverReportDigest = options.LinearSolverReportDigest;
report.linearSolverName = options.LinearSolverName;
report.linearSolverTolerance = options.LinearSolverTolerance;
report.linearSolverResidualNorm = options.LinearSolverResidualNorm;
report.linearSolverIterationCount = options.LinearSolverIterationCount;
report.expectedLinearSolverReportArtifactId = expectedLinearSolverReportArtifactId;
report.expectedLinearSolverReportDigest = expectedLinearSolverReportDigest;
report.expectedLinearSolverName = expectedLinearSolverName;
report.expectedLinearSolverTolerance = expectedLinearSolverTolerance;
report.expectedLinearSolverResidualNormMax = options.ExpectedLinearSolverResidualNormMax;
report.resultArtifactId = options.ResultArtifactId;
report.expectedResultArtifactId = expectedResultArtifactId;
report.runStartedAt = runStartedAt;
report.matlabVersion = matlabVersion;
report.expectedMatlabVersion = expectedMatlabVersion;
report.notebookSourceArtifactId = options.NotebookSourceArtifactId;
report.notebookSourceDigest = options.NotebookSourceDigest;
report.notebookSourcePath = options.NotebookSourcePath;
report.expectedNotebookSourceArtifactId = expectedNotebookSourceArtifactId;
report.expectedNotebookSourceDigest = expectedNotebookSourceDigest;
report.expectedNotebookSourcePath = expectedNotebookSourcePath;
report.parameterSetArtifactId = options.ParameterSetArtifactId;
report.parameterSetDigest = options.ParameterSetDigest;
report.parameterSetPath = options.ParameterSetPath;
report.expectedParameterSetArtifactId = expectedParameterSetArtifactId;
report.expectedParameterSetDigest = expectedParameterSetDigest;
report.expectedParameterSetPath = expectedParameterSetPath;
report.objectiveObservableId = options.ObjectiveObservableId;
report.objectiveObservableFamily = options.ObjectiveObservableFamily;
report.expectedObjectiveObservableId = expectedObjectiveObservableId;
report.expectedObjectiveObservableFamily = expectedObjectiveObservableFamily;
report.timingBreakdown = options.TimingBreakdown;
report.timingBreakdownNames = timingNames;
report.timingBreakdownSeconds = timingValues;
report.timingTotalSeconds = sum(timingValues);
report.execution = struct( ...
    "resultArtifactId", report.resultArtifactId, ...
    "linearSolverReportArtifactId", report.linearSolverReportArtifactId, ...
    "linearSolverReportDigest", report.linearSolverReportDigest, ...
    "linearSolverName", report.linearSolverName, ...
    "linearSolverTolerance", report.linearSolverTolerance, ...
    "linearSolverResidualNorm", report.linearSolverResidualNorm, ...
    "linearSolverIterationCount", report.linearSolverIterationCount, ...
    "runStartedAt", report.runStartedAt, ...
    "matlabVersion", report.matlabVersion, ...
    "notebookSourceArtifactId", report.notebookSourceArtifactId, ...
    "notebookSourceDigest", report.notebookSourceDigest, ...
    "notebookSourcePath", report.notebookSourcePath, ...
    "parameterSetArtifactId", report.parameterSetArtifactId, ...
    "parameterSetDigest", report.parameterSetDigest, ...
    "parameterSetPath", report.parameterSetPath, ...
    "objectiveObservableId", report.objectiveObservableId, ...
    "objectiveObservableFamily", report.objectiveObservableFamily, ...
    "timingBreakdown", report.timingBreakdown);
report.optimization = struct( ...
    "parameterSetArtifactId", report.parameterSetArtifactId, ...
    "parameterSetDigest", report.parameterSetDigest, ...
    "parameterSetPath", report.parameterSetPath, ...
    "objectiveObservableId", report.objectiveObservableId, ...
    "objectiveObservableFamily", report.objectiveObservableFamily);
report.requireTraceOutputArtifact = options.RequireTraceOutputArtifact;
report.requireNormalFluxArtifact = options.RequireNormalFluxArtifact;
report.requireLinearSolverReport = options.RequireLinearSolverReport;
report.requireNotebookSourceArtifact = options.RequireNotebookSourceArtifact;
report.requireParameterSetArtifact = options.RequireParameterSetArtifact;
report.requireCouplingConventionSchema = options.RequireCouplingConventionSchema;
report.requirePostprocessRowConventionSchema = options.RequirePostprocessRowConventionSchema;
report.requireTraceBasisSchema = options.RequireTraceBasisSchema;
report.operatorTraceSourceFileId = operatorTraceSourceFileId;
report.couplingKind = options.CouplingKind;
report.expectedCouplingKind = expectedCouplingKind;
report.formulationId = options.FormulationId;
report.expectedFormulationId = expectedFormulationId;
report.bemKernelFamily = options.BemKernelFamily;
report.expectedBemKernelFamily = expectedBemKernelFamily;
report.couplingConventionSchemaId = options.CouplingConventionSchemaId;
report.expectedCouplingConventionSchemaId = expectedCouplingConventionSchemaId;
report.postprocessRowConventionSchemaId = options.PostprocessRowConventionSchemaId;
report.expectedPostprocessRowConventionSchemaId = expectedPostprocessRowConventionSchemaId;
report.traceBasisSchemaId = traceBasisSchemaId;
report.operatorTraceBasisSchemaId = operatorTraceBasisSchemaId;
report.expectedTraceBasisSchemaId = expectedTraceBasisSchemaId;
report.bemKernelManifestId = options.BemKernelManifestId;
report.expectedBemKernelManifestId = expectedBemKernelManifestId;
report.bemKernelStrategy = options.BemKernelStrategy;
report.expectedBemKernelStrategy = expectedBemKernelStrategy;
report.kernelTimeConvention = options.KernelTimeConvention;
report.expectedKernelTimeConvention = expectedKernelTimeConvention;
report.assemblyRuleId = traceAssemblyRuleId;
report.operatorTraceAssemblyRuleId = operatorTraceAssemblyRuleId;
report.expectedAssemblyRuleId = expectedAssemblyRuleId;
report.quadratureRuleId = traceQuadratureRuleId;
report.operatorTraceQuadratureRuleId = operatorTraceQuadratureRuleId;
report.expectedQuadratureRuleId = expectedQuadratureRuleId;
report.volumeSpace = options.VolumeSpace;
report.expectedVolumeSpace = expectedVolumeSpace;
report.surfaceSpace = options.SurfaceSpace;
report.expectedSurfaceSpace = expectedSurfaceSpace;
report.boundaryNumbers = boundaryNumbers;
report.boundaryNames = boundaryNames;
report.boundaryRowIdentity = boundaryRowIdentity;
report.operatorBoundaryRowIdentity = operatorBoundaryRowIdentity;
report.expectedBoundaryNumbers = expectedBoundaryNumbers;
report.expectedBoundaryNames = expectedBoundaryNames;
report.polynomialOrder = 1;
report.curvedElementCount = 0;
report.geo = struct( ...
    "N", size(model.lukas.geo.nodes, 1), ...
    "conn_matrix", model.lukas.geo.conn_matrix);
report.gypsilab = struct( ...
    "elt", surfaceTriangles);
report.trace = struct( ...
    "fem_node_ids", traceNodeIds, ...
    "bem_node_ids", traceNodeIds, ...
    "source_file_id", traceSourceFileId, ...
    "boundary_numbers", boundaryNumbers, ...
    "boundary_names", boundaryNames, ...
    "boundary_row_identity", boundaryRowIdentity, ...
    "trace_row_identity", traceRowIdentity, ...
    "trace_matrix", T, ...
    "surface_triangles", surfaceTriangles, ...
    "surface_mesh_id", model.trace.surfaceMeshId, ...
    "trace_artifact_id", model.trace.traceArtifactId, ...
    "trace_operator_artifact_id", traceOperatorArtifactId, ...
    "trace_operator_policy", traceOperatorPolicy, ...
    "trace_output_artifact_id", traceOutputArtifactId, ...
    "trace_output_digest", traceOutputDigest, ...
    "trace_output_path", traceOutputPath, ...
    "trace_observable_id", traceObservableId, ...
    "trace_observable_family", traceObservableFamily, ...
    "normal_flux_artifact_id", options.NormalFluxArtifactId, ...
    "normal_flux_digest", options.NormalFluxDigest, ...
    "normal_flux_convention", options.NormalFluxConvention, ...
    "linear_solver_report_artifact_id", options.LinearSolverReportArtifactId, ...
    "linear_solver_report_digest", options.LinearSolverReportDigest, ...
    "linear_solver_name", options.LinearSolverName, ...
    "linear_solver_tolerance", options.LinearSolverTolerance, ...
    "linear_solver_residual_norm", options.LinearSolverResidualNorm, ...
    "linear_solver_iteration_count", options.LinearSolverIterationCount, ...
    "coupling_kind", options.CouplingKind, ...
    "formulation_id", options.FormulationId, ...
    "bem_kernel_family", options.BemKernelFamily, ...
    "coupling_convention_schema_id", options.CouplingConventionSchemaId, ...
    "fem_bem_postprocess_row_convention_schema_id", options.PostprocessRowConventionSchemaId, ...
    "trace_basis_schema_id", traceBasisSchemaId, ...
    "bem_kernel_manifest_id", options.BemKernelManifestId, ...
    "bem_kernel_strategy", options.BemKernelStrategy, ...
    "kernel_time_convention", options.KernelTimeConvention, ...
    "assembly_rule_id", traceAssemblyRuleId, ...
    "quadrature_rule_id", traceQuadratureRuleId, ...
    "volume_space", options.VolumeSpace, ...
    "surface_space", options.SurfaceSpace);

report.traceShape = size(T);
report.boundaryNodeIds = traceNodeIds;
report.traceRowIdentity = traceRowIdentity;
report.operatorTraceRowIdentity = operatorTraceRowIdentity;
report.boundaryTriangleCount = size(surfaceTriangles, 1);
report.volumeTetCount = size(model.lukas.geo.conn_matrix, 1);

report.checks = struct();
report.checks.meshIdRecorded = strlength(report.meshId) > 0;
report.checks.surfaceMeshIdRecorded = strlength(report.surfaceMeshId) > 0;
report.checks.traceArtifactIdRecorded = strlength(report.traceArtifactId) > 0;
report.checks.traceOperatorArtifactIdRecorded = strlength(report.traceOperatorArtifactId) > 0;
report.checks.operatorTraceOperatorArtifactIdRecorded = strlength(report.operatorTraceOperatorArtifactId) > 0;
report.checks.operatorTraceOperatorArtifactIdMatchesTrace = ...
    report.operatorTraceOperatorArtifactId == report.traceOperatorArtifactId;
report.checks.traceOperatorPolicyRecorded = strlength(report.traceOperatorPolicy) > 0;
report.checks.operatorTraceOperatorPolicyRecorded = strlength(report.operatorTraceOperatorPolicy) > 0;
report.checks.operatorTraceOperatorPolicyMatchesTrace = ...
    report.operatorTraceOperatorPolicy == report.traceOperatorPolicy;
report.checks.traceOperatorPolicyMatchesExpected = ...
    report.traceOperatorPolicy == expectedTraceOperatorPolicy;
report.checks.traceOperatorArtifactIdMatchesExpected = ...
    strlength(expectedTraceOperatorArtifactId) == 0 || ...
    report.traceOperatorArtifactId == expectedTraceOperatorArtifactId;
report.checks.traceOutputArtifactIdRecordedWhenRequired = ...
    ~report.requireTraceOutputArtifact || strlength(report.traceOutputArtifactId) > 0;
report.checks.traceOutputDigestRecordedWhenRequired = ...
    ~report.requireTraceOutputArtifact || strlength(report.traceOutputDigest) > 0;
report.checks.traceOutputPathRecordedWhenRequired = ...
    ~report.requireTraceOutputArtifact || strlength(report.traceOutputPath) > 0;
report.checks.operatorTraceOutputArtifactIdMatchesTrace = ...
    report.operatorTraceOutputArtifactId == report.traceOutputArtifactId;
report.checks.operatorTraceOutputDigestMatchesTrace = ...
    report.operatorTraceOutputDigest == report.traceOutputDigest;
report.checks.operatorTraceOutputPathMatchesTrace = ...
    report.operatorTraceOutputPath == report.traceOutputPath;
report.checks.traceOutputArtifactIdMatchesExpected = ...
    strlength(report.expectedTraceOutputArtifactId) == 0 || ...
    report.traceOutputArtifactId == report.expectedTraceOutputArtifactId;
report.checks.traceOutputDigestMatchesExpected = ...
    strlength(report.expectedTraceOutputDigest) == 0 || ...
    report.traceOutputDigest == report.expectedTraceOutputDigest;
report.checks.traceObservableIdRecorded = strlength(report.traceObservableId) > 0;
report.checks.operatorTraceObservableIdRecorded = strlength(report.operatorTraceObservableId) > 0;
report.checks.operatorTraceObservableIdMatchesTrace = ...
    report.operatorTraceObservableId == report.traceObservableId;
report.checks.traceObservableIdMatchesExpected = ...
    strlength(report.expectedTraceObservableId) == 0 || ...
    report.traceObservableId == report.expectedTraceObservableId;
report.checks.traceObservableFamilyRecorded = strlength(report.traceObservableFamily) > 0;
report.checks.operatorTraceObservableFamilyRecorded = strlength(report.operatorTraceObservableFamily) > 0;
report.checks.operatorTraceObservableFamilyMatchesTrace = ...
    report.operatorTraceObservableFamily == report.traceObservableFamily;
report.checks.traceObservableFamilyMatchesExpected = ...
    report.traceObservableFamily == expectedTraceObservableFamily;
report.checks.normalFluxConventionRecorded = strlength(report.normalFluxConvention) > 0;
report.checks.normalFluxConventionMatchesExpected = ...
    report.normalFluxConvention == expectedNormalFluxConvention;
report.checks.normalFluxArtifactIdRecordedWhenRequired = ...
    ~report.requireNormalFluxArtifact || strlength(report.normalFluxArtifactId) > 0;
report.checks.normalFluxDigestRecordedWhenRequired = ...
    ~report.requireNormalFluxArtifact || strlength(report.normalFluxDigest) > 0;
report.checks.normalFluxArtifactIdMatchesExpected = ...
    strlength(report.expectedNormalFluxArtifactId) == 0 || ...
    report.normalFluxArtifactId == report.expectedNormalFluxArtifactId;
report.checks.normalFluxDigestMatchesExpected = ...
    strlength(report.expectedNormalFluxDigest) == 0 || ...
    report.normalFluxDigest == report.expectedNormalFluxDigest;
report.checks.linearSolverReportArtifactIdRecordedWhenRequired = ...
    ~report.requireLinearSolverReport || strlength(report.linearSolverReportArtifactId) > 0;
report.checks.linearSolverReportDigestRecordedWhenRequired = ...
    ~report.requireLinearSolverReport || strlength(report.linearSolverReportDigest) > 0;
report.checks.linearSolverNameRecordedWhenRequired = ...
    ~report.requireLinearSolverReport || strlength(report.linearSolverName) > 0;
report.checks.linearSolverToleranceRecordedWhenRequired = ...
    ~report.requireLinearSolverReport || ~isnan(report.linearSolverTolerance);
report.checks.linearSolverResidualNormRecordedWhenRequired = ...
    ~report.requireLinearSolverReport || ~isnan(report.linearSolverResidualNorm);
report.checks.linearSolverToleranceFinitePositiveWhenPresent = ...
    isnan(report.linearSolverTolerance) || ...
    (isfinite(report.linearSolverTolerance) && report.linearSolverTolerance > 0);
report.checks.linearSolverResidualNormFiniteNonnegativeWhenPresent = ...
    isnan(report.linearSolverResidualNorm) || ...
    (isfinite(report.linearSolverResidualNorm) && report.linearSolverResidualNorm >= 0);
report.checks.linearSolverIterationCountNonnegativeWhenPresent = ...
    isnan(report.linearSolverIterationCount) || ...
    (isfinite(report.linearSolverIterationCount) && report.linearSolverIterationCount >= 0 && ...
    report.linearSolverIterationCount == fix(report.linearSolverIterationCount));
report.checks.linearSolverReportArtifactIdMatchesExpected = ...
    strlength(report.expectedLinearSolverReportArtifactId) == 0 || ...
    report.linearSolverReportArtifactId == report.expectedLinearSolverReportArtifactId;
report.checks.linearSolverReportDigestMatchesExpected = ...
    strlength(report.expectedLinearSolverReportDigest) == 0 || ...
    report.linearSolverReportDigest == report.expectedLinearSolverReportDigest;
report.checks.linearSolverNameMatchesExpected = ...
    strlength(report.expectedLinearSolverName) == 0 || ...
    report.linearSolverName == report.expectedLinearSolverName;
report.checks.linearSolverToleranceMatchesExpected = ...
    isnan(report.expectedLinearSolverTolerance) || ...
    (~isnan(report.linearSolverTolerance) && ...
    abs(report.linearSolverTolerance - report.expectedLinearSolverTolerance) <= ...
    max(1e-14, abs(report.expectedLinearSolverTolerance) * 1e-12));
report.checks.linearSolverResidualNormBelowExpectedMax = ...
    isnan(report.expectedLinearSolverResidualNormMax) || ...
    (~isnan(report.linearSolverResidualNorm) && ...
    report.linearSolverResidualNorm <= report.expectedLinearSolverResidualNormMax);
report.checks.resultArtifactIdRecorded = strlength(report.resultArtifactId) > 0;
report.checks.resultArtifactIdMatchesExpected = ...
    report.resultArtifactId == expectedResultArtifactId;
report.checks.runStartedAtRecorded = strlength(report.runStartedAt) > 0;
report.checks.runStartedAtIsoLike = isIsoTimestamp(report.runStartedAt);
report.checks.matlabVersionRecorded = strlength(report.matlabVersion) > 0;
report.checks.matlabVersionMatchesExpected = ...
    strlength(report.expectedMatlabVersion) == 0 || ...
    report.matlabVersion == report.expectedMatlabVersion;
report.checks.notebookSourceArtifactIdRecordedWhenRequired = ...
    ~report.requireNotebookSourceArtifact || strlength(report.notebookSourceArtifactId) > 0;
report.checks.notebookSourceDigestRecordedWhenRequired = ...
    ~report.requireNotebookSourceArtifact || strlength(report.notebookSourceDigest) > 0;
report.checks.notebookSourcePathRecordedWhenRequired = ...
    ~report.requireNotebookSourceArtifact || strlength(report.notebookSourcePath) > 0;
report.checks.notebookSourceArtifactIdMatchesExpected = ...
    strlength(report.expectedNotebookSourceArtifactId) == 0 || ...
    report.notebookSourceArtifactId == report.expectedNotebookSourceArtifactId;
report.checks.notebookSourceDigestMatchesExpected = ...
    strlength(report.expectedNotebookSourceDigest) == 0 || ...
    report.notebookSourceDigest == report.expectedNotebookSourceDigest;
report.checks.notebookSourcePathMatchesExpected = ...
    strlength(report.expectedNotebookSourcePath) == 0 || ...
    report.notebookSourcePath == report.expectedNotebookSourcePath;
report.checks.parameterSetArtifactIdRecordedWhenRequired = ...
    ~report.requireParameterSetArtifact || strlength(report.parameterSetArtifactId) > 0;
report.checks.parameterSetDigestRecordedWhenRequired = ...
    ~report.requireParameterSetArtifact || strlength(report.parameterSetDigest) > 0;
report.checks.parameterSetPathRecordedWhenRequired = ...
    ~report.requireParameterSetArtifact || strlength(report.parameterSetPath) > 0;
report.checks.parameterSetArtifactIdMatchesExpected = ...
    strlength(report.expectedParameterSetArtifactId) == 0 || ...
    report.parameterSetArtifactId == report.expectedParameterSetArtifactId;
report.checks.parameterSetDigestMatchesExpected = ...
    strlength(report.expectedParameterSetDigest) == 0 || ...
    report.parameterSetDigest == report.expectedParameterSetDigest;
report.checks.parameterSetPathMatchesExpected = ...
    strlength(report.expectedParameterSetPath) == 0 || ...
    report.parameterSetPath == report.expectedParameterSetPath;
report.checks.objectiveObservableIdMatchesExpected = ...
    strlength(report.expectedObjectiveObservableId) == 0 || ...
    report.objectiveObservableId == report.expectedObjectiveObservableId;
report.checks.objectiveObservableFamilyMatchesExpected = ...
    strlength(report.expectedObjectiveObservableFamily) == 0 || ...
    report.objectiveObservableFamily == report.expectedObjectiveObservableFamily;
report.checks.timingBreakdownRecorded = numel(report.timingBreakdownNames) > 0;
report.checks.timingBreakdownHasFourItems = numel(report.timingBreakdownNames) >= 4;
report.checks.timingBreakdownHasAtMostFourItems = numel(report.timingBreakdownNames) <= 4;
report.checks.timingBreakdownFiniteNonnegative = ...
    all(isfinite(report.timingBreakdownSeconds)) && ...
    all(report.timingBreakdownSeconds >= 0);
report.checks.sourceFileIdRecorded = strlength(report.sourceFileId) > 0;
report.checks.traceSourceFileIdRecorded = strlength(string(report.trace.source_file_id)) > 0;
report.checks.traceSourceFileIdMatchesIdentity = string(report.trace.source_file_id) == report.sourceFileId;
report.checks.operatorTraceSourceFileIdRecorded = strlength(report.operatorTraceSourceFileId) > 0;
report.checks.operatorTraceSourceFileIdMatchesIdentity = report.operatorTraceSourceFileId == report.sourceFileId;
report.checks.boundaryNumbersRecorded = numel(boundaryNumbers) == size(surfaceTriangles, 1);
report.checks.boundaryNamesRecorded = numel(boundaryNames) == size(surfaceTriangles, 1) && all(strlength(boundaryNames) > 0);
report.checks.operatorTraceBoundaryNumbersMatch = isequal(operatorBoundaryNumbers(:), boundaryNumbers(:));
report.checks.operatorTraceBoundaryNamesMatch = isequal(operatorBoundaryNames(:), boundaryNames(:));
report.checks.boundaryRowIdentityRecorded = numel(boundaryRowIdentity) == size(surfaceTriangles, 1);
report.checks.boundaryRowIdentityRowIndicesMatch = isequal(boundaryRowIndices, (1:size(surfaceTriangles, 1)).');
report.checks.boundaryRowIdentityTrianglesMatch = isequal(boundaryRowTriangles, surfaceTriangles);
report.checks.boundaryRowIdentityNumbersMatch = isequal(boundaryRowNumbers(:), boundaryNumbers(:));
report.checks.boundaryRowIdentityNamesMatch = isequal(boundaryRowNames(:), boundaryNames(:));
report.checks.operatorTraceBoundaryRowIdentityRecorded = numel(operatorBoundaryRowIdentity) == size(surfaceTriangles, 1);
report.checks.operatorTraceBoundaryRowIdentityMatchesTrace = isequal(operatorBoundaryRowIdentity, boundaryRowIdentity);
report.checks.operatorTraceBoundaryRowIdentityRowIndicesMatch = isequal( ...
    operatorBoundaryRowIndices, (1:size(surfaceTriangles, 1)).');
report.checks.operatorTraceBoundaryRowIdentityTrianglesMatch = isequal(operatorBoundaryRowTriangles, surfaceTriangles);
report.checks.operatorTraceBoundaryRowIdentityNumbersMatch = isequal(operatorBoundaryRowNumbers(:), boundaryNumbers(:));
report.checks.operatorTraceBoundaryRowIdentityNamesMatch = isequal(operatorBoundaryRowNames(:), boundaryNames(:));
report.checks.boundaryNumbersMatchExpected = isempty(expectedBoundaryNumbers) || ...
    isequal(sort(unique(boundaryNumbers(:))), expectedBoundaryNumbers);
report.checks.boundaryNamesMatchExpected = isempty(expectedBoundaryNames) || ...
    isequal(sort(unique(boundaryNames(:))), expectedBoundaryNames);
report.checks.sourceFormatIsVol = report.sourceFormat == ".vol";
report.checks.volumeSpaceMatchesExpected = report.volumeSpace == expectedVolumeSpace;
report.checks.surfaceSpaceMatchesExpected = report.surfaceSpace == expectedSurfaceSpace;
report.checks.couplingKindMatchesExpected = report.couplingKind == expectedCouplingKind;
report.checks.formulationIdMatchesExpected = report.formulationId == expectedFormulationId;
report.checks.bemKernelFamilyMatchesExpected = report.bemKernelFamily == expectedBemKernelFamily;
report.checks.couplingConventionSchemaIdRecordedWhenRequired = ...
    ~report.requireCouplingConventionSchema || strlength(report.couplingConventionSchemaId) > 0;
report.checks.couplingConventionSchemaIdMatchesExpected = ...
    strlength(report.expectedCouplingConventionSchemaId) == 0 || ...
    report.couplingConventionSchemaId == report.expectedCouplingConventionSchemaId;
report.checks.postprocessRowConventionSchemaIdRecordedWhenRequired = ...
    ~report.requirePostprocessRowConventionSchema || strlength(report.postprocessRowConventionSchemaId) > 0;
report.checks.postprocessRowConventionSchemaIdMatchesExpected = ...
    strlength(report.expectedPostprocessRowConventionSchemaId) == 0 || ...
    report.postprocessRowConventionSchemaId == report.expectedPostprocessRowConventionSchemaId;
report.checks.traceBasisSchemaIdRecordedWhenRequired = ...
    ~report.requireTraceBasisSchema || strlength(report.traceBasisSchemaId) > 0;
report.checks.traceBasisSchemaIdMatchesExpected = ...
    strlength(report.expectedTraceBasisSchemaId) == 0 || ...
    report.traceBasisSchemaId == report.expectedTraceBasisSchemaId;
report.checks.operatorTraceBasisSchemaIdRecorded = strlength(report.operatorTraceBasisSchemaId) > 0;
report.checks.operatorTraceBasisSchemaIdMatchesTrace = ...
    report.operatorTraceBasisSchemaId == report.traceBasisSchemaId;
report.checks.bemKernelManifestIdRecorded = strlength(report.bemKernelManifestId) > 0;
report.checks.bemKernelManifestIdMatchesExpected = ...
    report.bemKernelManifestId == expectedBemKernelManifestId;
report.checks.bemKernelStrategyRecorded = strlength(report.bemKernelStrategy) > 0;
report.checks.bemKernelStrategyMatchesExpected = ...
    report.bemKernelStrategy == expectedBemKernelStrategy;
report.checks.kernelTimeConventionRecorded = strlength(report.kernelTimeConvention) > 0;
report.checks.kernelTimeConventionMatchesExpected = ...
    report.kernelTimeConvention == expectedKernelTimeConvention;
report.checks.assemblyRuleIdRecorded = strlength(report.assemblyRuleId) > 0;
report.checks.operatorTraceAssemblyRuleIdRecorded = strlength(report.operatorTraceAssemblyRuleId) > 0;
report.checks.operatorTraceAssemblyRuleIdMatchesTrace = ...
    report.operatorTraceAssemblyRuleId == report.assemblyRuleId;
report.checks.assemblyRuleIdMatchesExpected = ...
    report.assemblyRuleId == expectedAssemblyRuleId;
report.checks.quadratureRuleIdRecorded = strlength(report.quadratureRuleId) > 0;
report.checks.operatorTraceQuadratureRuleIdRecorded = strlength(report.operatorTraceQuadratureRuleId) > 0;
report.checks.operatorTraceQuadratureRuleIdMatchesTrace = ...
    report.operatorTraceQuadratureRuleId == report.quadratureRuleId;
report.checks.quadratureRuleIdMatchesExpected = ...
    report.quadratureRuleId == expectedQuadratureRuleId;
report.checks.traceRowsMatchBoundaryNodes = size(T, 1) == numel(traceNodeIds);
report.checks.traceColumnsMatchVolumeNodes = size(T, 2) == size(model.lukas.geo.nodes, 1);
report.checks.traceRowsAreOneHot = all(rowNonzeroCount == 1) && all(abs(rowMax - 1) < 1e-14);
report.checks.traceRowIdentityRecorded = numel(traceRowIdentity) == numel(traceNodeIds);
report.checks.traceRowIdentityRowIndicesMatch = isequal(traceRowIndex, (1:numel(traceNodeIds)).');
report.checks.traceRowIdentityFemNodesMatch = isequal(traceRowFemNodeIds, traceNodeIds);
report.checks.traceRowIdentityBemNodesMatch = isequal(traceRowBemNodeIds, traceNodeIds);
report.checks.traceRowIdentityUnique = ...
    numel(unique(traceRowFemNodeIds)) == numel(traceRowFemNodeIds) && ...
    numel(unique(traceRowBemNodeIds)) == numel(traceRowBemNodeIds);
report.checks.traceRowIdentityMatchesTraceMatrix = isequal(traceRowFemNodeIds, traceMatrixNodeIds);
report.checks.operatorTraceRowIdentityRecorded = numel(operatorTraceRowIdentity) == numel(traceNodeIds);
report.checks.operatorTraceRowIdentityMatchesTrace = isequal(operatorTraceRowIdentity, traceRowIdentity);
report.checks.operatorTraceRowIdentityRowIndicesMatch = isequal(operatorTraceRowIndex, (1:numel(traceNodeIds)).');
report.checks.operatorTraceRowIdentityFemNodesMatch = isequal(operatorTraceRowFemNodeIds, traceNodeIds);
report.checks.operatorTraceRowIdentityBemNodesMatch = isequal(operatorTraceRowBemNodeIds, traceNodeIds);
report.checks.operatorTraceRowIdentityMatchesTraceMatrix = isequal(operatorTraceRowFemNodeIds, traceMatrixNodeIds);

if all(structfun(@(value) logical(value), report.checks))
    report.status = "ok";
else
    report.status = "needs_attention";
end

report.notes = [
    "Use this manifest before comparing FEM/BEM coupled values."
    "Trace, surface mesh, formulation, kernel family, and spaces must be one package."
    "The FEM/BEM coupling convention schema id binds normals, trace rows, spaces, kernel family, and formulation into one reusable teaching contract."
    "The postprocess row convention schema id binds trace least-squares rows, residual/objective reduction, and optimizer-facing scalar rows separately from the coupling convention."
    "The trace basis schema id binds H1 P1 volume nodal basis rows to scalar P1 surface/BEM basis rows before values are reused."
    "The BEM kernel manifest id, strategy, and time convention must travel with the coupling package before low-frequency acoustic rows are reused."
    "The assembly rule id and quadrature rule id must travel with the operator rows before notebook results are reused."
    "The .vol source file identity must match the trace artifact identity before reuse."
    "The trace operator artifact must be separate from later observation or field-map operators."
    "The trace output artifact id, digest, and path identify the concrete matrix consumed by notebooks or later coupling steps."
    "The trace observable id and family identify this result as a FEM/BEM boundary trace rather than a remote field map, residual curve, or optimizer diagnostic."
    "Normal-flux artifact id, digest, and convention identify the orientation/sign evidence consumed by Neumann, acoustic, or surface-source rows."
    "Linear solver report artifact id, digest, solver name, tolerance, residual norm, and iteration count identify the concrete solve used by the result package."
    "Result artifact id, MATLAB version, ISO-like run timestamp, and a compact four-stage timing breakdown should travel with executed notebooks and JSON results."
    "Notebook/source artifact id, digest, and path identify the teaching code that produced the executed result package."
    "Parameter-set artifact id, digest, and path identify the initial values or design variables reused as notebook defaults or optimization inputs."
    "Objective observable id and family identify the scalar quantity that an optimizer, panel, or replayed notebook will consume."
    "Boundary number and boundary name identity must match the trace before boundary-condition rows reuse the package."
    "Boundary row identity binds each surface triangle row to its nodes, boundary number, and boundary name before BEM-kernel rows reuse it."
    "The trace operator records its row identity explicitly so notebooks do not need to infer it from sparse nonzeros."
    "This is a readable first-order teaching contract, not a high-performance solver API."
    ];
end


function tf = isIsoTimestamp(value)
text = char(string(value));
tf = ~isempty(regexp(text, ...
    '^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(Z|[+-]\d{2}:\d{2})?$', ...
    'once'));
end


function identity = makeTraceRowIdentity(traceNodeIds)
nRows = numel(traceNodeIds);
if nRows == 0
    identity = struct( ...
        "trace_row_index", {}, ...
        "fem_node_id", {}, ...
        "bem_node_id", {}, ...
        "surface_node_index", {});
    return
end

identity = struct( ...
    "trace_row_index", num2cell((1:nRows).'), ...
    "fem_node_id", num2cell(traceNodeIds(:)), ...
    "bem_node_id", num2cell(traceNodeIds(:)), ...
    "surface_node_index", num2cell((1:nRows).'));
end


function expected = defaultExpected(value, fallback)
if strlength(value) == 0
    expected = fallback;
else
    expected = value;
end
end


function expected = defaultExpectedDouble(value, fallback)
if isnan(value)
    expected = fallback;
else
    expected = value;
end
end


function value = getStringField(record, name)
if isstruct(record) && isfield(record, name)
    value = string(record.(name));
else
    value = "";
end
end


function value = getNumericVectorField(record, name)
if isstruct(record) && isfield(record, name)
    value = record.(name);
    value = value(:);
else
    value = [];
end
end


function value = getStringVectorField(record, name)
if isstruct(record) && isfield(record, name)
    value = string(record.(name));
    value = value(:);
else
    value = strings(0, 1);
end
end


function value = getStructVectorField(record, name)
if isstruct(record) && isfield(record, name)
    value = record.(name);
    value = value(:);
else
    value = struct( ...
        "trace_row_index", {}, ...
        "fem_node_id", {}, ...
        "bem_node_id", {}, ...
        "surface_node_index", {});
end
end


function value = getBoundaryRowIdentityField(record, name)
if isstruct(record) && isfield(record, name)
    value = record.(name);
    value = value(:);
else
    value = struct( ...
        "surface_triangle_index", {}, ...
        "surface_triangle_nodes", {}, ...
        "boundary_number", {}, ...
        "boundary_name", {}, ...
        "adjacent_tet_index", {});
end
end


function [names, values] = timingBreakdownArrays(timing)
names = string(fieldnames(timing));
values = zeros(numel(names), 1);
for k = 1:numel(names)
    values(k) = double(timing.(names(k)));
end
end


function [indices, triangles, numbers, names] = boundaryRowIdentityArrays(identity)
nRows = numel(identity);
indices = zeros(nRows, 1);
triangles = zeros(nRows, 3);
numbers = zeros(nRows, 1);
names = strings(nRows, 1);
for k = 1:nRows
    indices(k) = identity(k).surface_triangle_index;
    triangles(k, :) = identity(k).surface_triangle_nodes;
    numbers(k) = identity(k).boundary_number;
    names(k) = string(identity(k).boundary_name);
end
end
