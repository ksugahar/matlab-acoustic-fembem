function report = femBemCouplingManifest(model, options)
%FEMBEMCOUPLINGMANIFEST Record readable FEM/BEM coupling identity.
%
% This helper stops before solving.  It makes the trace handoff explicit so
% students can see which volume space, surface space, formulation, and BEM
% kernel family belong to the same .vol-derived FEM/BEM package.
%
% The manifest RECORDS the package; femBemManifestChecks VALIDATES the
% record.  Boundary-condition identity is read once from the SurfaceMesh;
% the trace identity is read from both the model TraceOperator and the
% assembled operators bundle so a stale operator bundle is caught.

arguments
    model (1,1) FemBemModel
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

if isempty(model.operators)
    model = model.assemble();
end
trace = model.trace;               % TraceOperator held by the model
operatorTrace = model.operators.trace;   % TraceOperator inside the bundle

% --- resolve expected values (empty expected = "expect the recorded value") ---
expectedCouplingKind = defaultExpected(options.ExpectedCouplingKind, options.CouplingKind);
expectedFormulationId = defaultExpected(options.ExpectedFormulationId, options.FormulationId);
expectedBemKernelFamily = defaultExpected(options.ExpectedBemKernelFamily, options.BemKernelFamily);
expectedCouplingConventionSchemaId = defaultExpected( ...
    options.ExpectedCouplingConventionSchemaId, ...
    options.CouplingConventionSchemaId);
expectedPostprocessRowConventionSchemaId = defaultExpected( ...
    options.ExpectedPostprocessRowConventionSchemaId, ...
    options.PostprocessRowConventionSchemaId);
traceBasisSchemaId = defaultExpected(options.TraceBasisSchemaId, trace.basisSchemaId);
expectedTraceBasisSchemaId = defaultExpected( ...
    options.ExpectedTraceBasisSchemaId, ...
    traceBasisSchemaId);
expectedBemKernelManifestId = defaultExpected( ...
    options.ExpectedBemKernelManifestId, options.BemKernelManifestId);
expectedBemKernelStrategy = defaultExpected( ...
    options.ExpectedBemKernelStrategy, options.BemKernelStrategy);
expectedKernelTimeConvention = defaultExpected( ...
    options.ExpectedKernelTimeConvention, options.KernelTimeConvention);
expectedAssemblyRuleId = defaultExpected( ...
    options.ExpectedAssemblyRuleId, options.AssemblyRuleId);
expectedQuadratureRuleId = defaultExpected( ...
    options.ExpectedQuadratureRuleId, options.QuadratureRuleId);
expectedVolumeSpace = defaultExpected(options.ExpectedVolumeSpace, options.VolumeSpace);
expectedSurfaceSpace = defaultExpected(options.ExpectedSurfaceSpace, options.SurfaceSpace);
expectedTraceOperatorPolicy = defaultExpected( ...
    options.ExpectedTraceOperatorPolicy, ...
    "one_hot_boundary_node_injection_from_vol_node_ids");
expectedTraceObservableFamily = defaultExpected( ...
    options.ExpectedTraceObservableFamily, ...
    "fem_bem_boundary_trace");
expectedNormalFluxConvention = defaultExpected( ...
    options.ExpectedNormalFluxConvention, options.NormalFluxConvention);
expectedLinearSolverReportArtifactId = defaultExpected( ...
    options.ExpectedLinearSolverReportArtifactId, options.LinearSolverReportArtifactId);
expectedLinearSolverReportDigest = defaultExpected( ...
    options.ExpectedLinearSolverReportDigest, options.LinearSolverReportDigest);
expectedLinearSolverName = defaultExpected( ...
    options.ExpectedLinearSolverName, options.LinearSolverName);
expectedLinearSolverTolerance = defaultExpectedDouble( ...
    options.ExpectedLinearSolverTolerance, options.LinearSolverTolerance);
expectedResultArtifactId = defaultExpected( ...
    options.ExpectedResultArtifactId, options.ResultArtifactId);
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

% --- read the coupled package off the model objects ---
T = full(trace.matrix);
traceNodeIds = trace.femNodeIds(:);
bemNodeIds = trace.bemNodeIds(:);
surfaceTriangles = model.surface.triGlobal;
boundaryNumbers = model.surface.col(:);
boundaryNames = model.surface.names(:);
boundaryRowIdentity = model.surface.rowIdentity(:);
traceRowIdentity = makeTraceRowIdentity(traceNodeIds);
[~, ~, sourceExt] = fileparts(model.mesh.sourcePath);

report = struct();
report.policy = "readable_fem_bem_coupling_manifest_gate";
report.meshId = model.mesh.meshId;
report.exportId = model.mesh.meshId;
report.sourcePath = model.mesh.sourcePath;
report.sourceFileId = model.mesh.sourceFileId;
report.sourceFormat = string(sourceExt);
report.surfaceMeshId = trace.surfaceMeshId;
report.traceArtifactId = trace.artifactId;
report.traceOperatorArtifactId = trace.operatorArtifactId;
report.operatorTraceOperatorArtifactId = operatorTrace.operatorArtifactId;
report.traceOperatorPolicy = trace.operatorPolicy;
report.operatorTraceOperatorPolicy = operatorTrace.operatorPolicy;
report.expectedTraceOperatorArtifactId = options.ExpectedTraceOperatorArtifactId;
report.expectedTraceOperatorPolicy = expectedTraceOperatorPolicy;
report.traceOutputArtifactId = trace.outputArtifactId;
report.operatorTraceOutputArtifactId = operatorTrace.outputArtifactId;
report.traceOutputDigest = trace.outputDigest;
report.operatorTraceOutputDigest = operatorTrace.outputDigest;
report.traceOutputPath = trace.outputPath;
report.operatorTraceOutputPath = operatorTrace.outputPath;
report.expectedTraceOutputArtifactId = options.ExpectedTraceOutputArtifactId;
report.expectedTraceOutputDigest = options.ExpectedTraceOutputDigest;
report.traceObservableId = trace.observableId;
report.operatorTraceObservableId = operatorTrace.observableId;
report.expectedTraceObservableId = options.ExpectedTraceObservableId;
report.traceObservableFamily = trace.observableFamily;
report.operatorTraceObservableFamily = operatorTrace.observableFamily;
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
report.expectedMatlabVersion = options.ExpectedMatlabVersion;
report.notebookSourceArtifactId = options.NotebookSourceArtifactId;
report.notebookSourceDigest = options.NotebookSourceDigest;
report.notebookSourcePath = options.NotebookSourcePath;
report.expectedNotebookSourceArtifactId = options.ExpectedNotebookSourceArtifactId;
report.expectedNotebookSourceDigest = options.ExpectedNotebookSourceDigest;
report.expectedNotebookSourcePath = options.ExpectedNotebookSourcePath;
report.parameterSetArtifactId = options.ParameterSetArtifactId;
report.parameterSetDigest = options.ParameterSetDigest;
report.parameterSetPath = options.ParameterSetPath;
report.expectedParameterSetArtifactId = options.ExpectedParameterSetArtifactId;
report.expectedParameterSetDigest = options.ExpectedParameterSetDigest;
report.expectedParameterSetPath = options.ExpectedParameterSetPath;
report.objectiveObservableId = options.ObjectiveObservableId;
report.objectiveObservableFamily = options.ObjectiveObservableFamily;
report.expectedObjectiveObservableId = options.ExpectedObjectiveObservableId;
report.expectedObjectiveObservableFamily = options.ExpectedObjectiveObservableFamily;
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
report.operatorTraceSourceFileId = operatorTrace.sourceFileId;
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
report.operatorTraceBasisSchemaId = operatorTrace.basisSchemaId;
report.expectedTraceBasisSchemaId = expectedTraceBasisSchemaId;
report.bemKernelManifestId = options.BemKernelManifestId;
report.expectedBemKernelManifestId = expectedBemKernelManifestId;
report.bemKernelStrategy = options.BemKernelStrategy;
report.expectedBemKernelStrategy = expectedBemKernelStrategy;
report.kernelTimeConvention = options.KernelTimeConvention;
report.expectedKernelTimeConvention = expectedKernelTimeConvention;
report.assemblyRuleId = trace.assemblyRuleId;
report.operatorTraceAssemblyRuleId = operatorTrace.assemblyRuleId;
report.expectedAssemblyRuleId = expectedAssemblyRuleId;
report.quadratureRuleId = trace.quadratureRuleId;
report.operatorTraceQuadratureRuleId = operatorTrace.quadratureRuleId;
report.expectedQuadratureRuleId = expectedQuadratureRuleId;
report.volumeSpace = options.VolumeSpace;
report.expectedVolumeSpace = expectedVolumeSpace;
report.surfaceSpace = options.SurfaceSpace;
report.expectedSurfaceSpace = expectedSurfaceSpace;
report.boundaryNumbers = boundaryNumbers;
report.boundaryNames = boundaryNames;
report.boundaryRowIdentity = boundaryRowIdentity;
report.expectedBoundaryNumbers = sort(unique(options.ExpectedBoundaryNumbers(:)));
report.expectedBoundaryNames = sort(unique(options.ExpectedBoundaryNames(:)));
report.polynomialOrder = 1;
report.curvedElementCount = 0;
report.geo = struct( ...
    "N", size(model.mesh.vtx, 1), ...
    "conn_matrix", model.mesh.tet);
report.gypsilab = struct( ...
    "elt", surfaceTriangles);
report.trace = struct( ...
    "fem_node_ids", traceNodeIds, ...
    "bem_node_ids", bemNodeIds, ...
    "source_file_id", trace.sourceFileId, ...
    "boundary_numbers", boundaryNumbers, ...
    "boundary_names", boundaryNames, ...
    "boundary_row_identity", boundaryRowIdentity, ...
    "trace_row_identity", traceRowIdentity, ...
    "trace_matrix", T, ...
    "surface_triangles", surfaceTriangles, ...
    "surface_mesh_id", trace.surfaceMeshId, ...
    "trace_artifact_id", trace.artifactId, ...
    "trace_operator_artifact_id", trace.operatorArtifactId, ...
    "trace_operator_policy", trace.operatorPolicy, ...
    "trace_output_artifact_id", trace.outputArtifactId, ...
    "trace_output_digest", trace.outputDigest, ...
    "trace_output_path", trace.outputPath, ...
    "trace_observable_id", trace.observableId, ...
    "trace_observable_family", trace.observableFamily, ...
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
    "assembly_rule_id", trace.assemblyRuleId, ...
    "quadrature_rule_id", trace.quadratureRuleId, ...
    "volume_space", options.VolumeSpace, ...
    "surface_space", options.SurfaceSpace);

report.traceShape = size(T);
report.boundaryNodeIds = traceNodeIds;
report.traceRowIdentity = traceRowIdentity;
report.operatorTraceRowIdentity = operatorTrace.rowIdentity(:);
report.boundaryTriangleCount = size(surfaceTriangles, 1);
report.volumeTetCount = size(model.mesh.tet, 1);

report.checks = femBemManifestChecks(report);

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
    "bem_node_id", num2cell((1:nRows).'), ...
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


function [names, values] = timingBreakdownArrays(timing)
names = string(fieldnames(timing));
values = zeros(numel(names), 1);
for k = 1:numel(names)
    values(k) = double(timing.(names(k)));
end
end
