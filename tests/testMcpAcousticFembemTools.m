function tests = testMcpAcousticFembemTools
%TESTMCPACOUSTICFEMBEMTOOLS MCP-facing acoustic FEM/BEM entry points.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(repoRoot));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));
end


function testRepositoryRoot(testCase)
root = acoustic_fembem.repository_root();
verifyTrue(testCase, isfolder(fullfile(root, "matlab_api")));
verifyTrue(testCase, isfile(fullfile(root, "mcp", "extensions", "acoustic-fembem-tools.json")));
end


function testMcpLayerStaysInsideAcousticFembemRepo(testCase)
root = acoustic_fembem.repository_root();
readme = string(fileread(fullfile(root, "mcp", "README.md")));
requirements = string(fileread(fullfile(root, "mcp", "REQUIREMENTS.md")));
verifySubstring(testCase, readme, "stays inside this repository");
verifySubstring(testCase, readme, "official MathWorks server remains the runtime");
verifySubstring(testCase, readme, "MATLAB Agentic Toolkit");
verifySubstring(testCase, readme, "existing-session workflow");
verifySubstring(testCase, requirements, "extension is intentionally thin");
verifySubstring(testCase, requirements, "process/session management");
verifySubstring(testCase, requirements, "matlab-agentic-toolkit");
verifySubstring(testCase, requirements, "Existing MATLAB Session Policy");
end


function testKnowledgeIncludesCrossvalTopic(testCase)
body = acoustic_fembem.fembem_knowledge("radia_ngsolve_crossval");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, ".vol");
verifySubstring(testCase, body, "radia-ngsolve");
end


function testKnowledgeIncludesNgsolveBem50Topic(testCase)
body = acoustic_fembem.fembem_knowledge("ngsolve_bem_50");
verifyGreaterThan(testCase, strlength(body), 600);
verifySubstring(testCase, body, "50 examples");
verifySubstring(testCase, body, "NGSolve.BEM");
verifySubstring(testCase, body, "Netgen");
verifySubstring(testCase, body, "Do not claim an analytic solution");
end


function testKnowledgeIncludesPdeVolBridgeTopic(testCase)
body = acoustic_fembem.fembem_knowledge("pde_vol_bridge");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, "PDE Toolbox");
verifySubstring(testCase, body, "writePdeMeshVol");
verifySubstring(testCase, body, ".vol");
end


function testKnowledgeIncludesVolVisualizationTopic(testCase)
body = acoustic_fembem.fembem_knowledge("vol_visualization");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, "Netgen");
verifySubstring(testCase, body, "plotVolMesh");
verifySubstring(testCase, body, "acoustic_fembem_vol_mesh_summary");
verifySubstring(testCase, body, "Ng_LoadMesh");
verifySubstring(testCase, body, ".sol files are mesh-free");
end


function testKnowledgeIncludesConvolutionQuadratureTopic(testCase)
body = acoustic_fembem.fembem_knowledge("convolution_quadrature");
verifyGreaterThan(testCase, strlength(body), 400);
verifySubstring(testCase, body, "Lubich CQ");
verifySubstring(testCase, body, "laplaceSingleLayerGalerkin");
verifySubstring(testCase, body, "A-stability");
verifySubstring(testCase, body, "IMAGINARY axis");
verifyEqual(testCase, acoustic_fembem.fembem_knowledge("cq"), body);
end


function testCqTimeGridMcpContract(testCase)
report = cqTimeGridManifest(1e-3, 100, "BDF2", 128, 0);
verifyTrue(testCase, report.ok);
verifyEqual(testCase, report.contourSamples, 128);
verifyEqual(testCase, report.timeEnd, 0.099, "AbsTol", 1e-14);
verifyGreaterThan(testCase, report.minRealLaplaceNode, 0);

bad = cqTimeGridManifest(1e-3, 100, "BDF2", 64, 0.9);
verifyFalse(testCase, bad.ok);
verifyFalse(testCase, bad.checks.contourCoversTimeGrid);

out = evalc("acoustic_fembem.check_cq_time_grid(1e-3, 100, ""BDF2"", 128, 0)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_cq_time_grid");
verifyEqual(testCase, string(decoded.result.schema), "matlab-acoustic-fembem.cq-time-grid.v1");
end


function testBalancedLearningProfileMcpContract(testCase)
profile = balancedLearningProfile();
verifyEqual(testCase, profile.policy, "equal_capability_gain_v1");
verifyEqual(testCase, profile.stage_count, 10);
verifyEqual(testCase, numel(unique([profile.stages.capability_id])), 10);
verifyEqual(testCase, sort(string(fieldnames(profile.workflow_roles))), ...
    sort(["detect"; "check"; "run"; "test"]));
verifyEqual(testCase, profile.self_check.status, "ok");

bad = profile;
bad.stages(4).negative_control = "";
rejected = validateBalancedLearningProfile(bad);
verifyEqual(testCase, rejected.status, "needs_attention");
verifyFalse(testCase, rejected.checks.controls_complete);

out = evalc("acoustic_fembem.check_balanced_learning_profile()");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_balanced_learning_profile");
end


function testKnowledgeIncludesPublicAcousticBlogLessons(testCase)
body = acoustic_fembem.fembem_knowledge("public_acoustic_blog_lessons");
verifyGreaterThan(testCase, strlength(body), 700);
verifySubstring(testCase, body, "Unbounded exterior radiation");
verifySubstring(testCase, body, "frequency domain");
verifySubstring(testCase, body, "high-order surface");
verifySubstring(testCase, body, "PML=false");
verifySubstring(testCase, body, "time-domain lane is explicitly CQ");
verifySubstring(testCase, body, "Acoustic-structure interaction");
verifySubstring(testCase, body, "two-way coupling");
verifySubstring(testCase, body, "Impedance lumping");
verifySubstring(testCase, body, "p=Z_s v");
verifySubstring(testCase, body, "p=Z Q");
verifySubstring(testCase, body, "local reaction and extended reaction");
verifySubstring(testCase, body, "high-order Zs");
verifySubstring(testCase, body, "Schroeder-frequency");
verifySubstring(testCase, body, "acoustic_method_selection_manifest_gate");
end


function testKnowledgeIncludesPublicAcousticNonboundary10(testCase)
body = acoustic_fembem.fembem_knowledge("public_acoustic_nonboundary_10");
verifyGreaterThan(testCase, strlength(body), 1200);
verifySubstring(testCase, body, "10 public acoustic non-boundary problems");
verifySubstring(testCase, body, "absorbing boundaries");
verifySubstring(testCase, body, "PML=false");
verifySubstring(testCase, body, "Acoustic trap");
verifySubstring(testCase, body, "Surface-acoustic-wave droplet streaming");
verifySubstring(testCase, body, "Thermoviscous acoustic radiation force");
verifySubstring(testCase, body, "Thermoacoustic engine");
verifySubstring(testCase, body, "Acoustic topology optimization");
verifySubstring(testCase, body, "Room response split");
verifySubstring(testCase, body, "Small-speaker room impulse response");
verifySubstring(testCase, body, "Ultrasonic pipe pulse-echo");
verifySubstring(testCase, body, "public_acoustic_nonboundary_problem_catalog");
verifySubstring(testCase, body, "acoustic_nonboundary_problem_catalog_manifest_gate");
end


function testKnowledgeIncludesGmshArtifactTopic(testCase)
% The MATLAB lane is gmsh-free: the topic documents the native GIF path and
% points gmsh acoustic movies to the radia-acoustic (GmshPostExport) side.
body = acoustic_fembem.fembem_knowledge("gmsh_artifact");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "gmsh-free");
verifySubstring(testCase, body, "writeSoftSphereScatterGif");
verifySubstring(testCase, body, "drumScatterField");
verifySubstring(testCase, body, "does not require Gmsh");
verifySubstring(testCase, body, "GmshPostExport");
verifySubstring(testCase, body, "radia-acoustic");
verifySubstring(testCase, body, ".msh v4.1");
end


function testKnowledgeIncludesCatalog100Topic(testCase)
body = acoustic_fembem.fembem_knowledge("catalog_100");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "100-case");
verifySubstring(testCase, body, "GYP-001..010");
verifySubstring(testCase, body, "GYP-091..100");
end


function testKnowledgeIncludesDrumTopic(testCase)
body = acoustic_fembem.fembem_knowledge("vibroacoustic_drum");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "baffled circular membrane");
verifySubstring(testCase, body, "normal velocity");
verifySubstring(testCase, body, "NGSolve.BEM");
verifySubstring(testCase, body, "the drum structure is FEM");
verifySubstring(testCase, body, "air radiation");
verifySubstring(testCase, body, "acoustic BEM");
verifySubstring(testCase, body, "plotDrumStepTimeField");
verifySubstring(testCase, body, "writeDrumStepTimeGif");
verifySubstring(testCase, body, "drumHighOrderImpedanceScene");
verifySubstring(testCase, body, "drumFemBemCoupledDemo");
verifySubstring(testCase, body, "volFemBemIfftResponse");
verifySubstring(testCase, body, "volTdBemConvolutionQuadrature");
verifySubstring(testCase, body, "axis-equal");
verifySubstring(testCase, body, "not a hemisphere");
verifySubstring(testCase, body, "top membrane");
verifySubstring(testCase, body, "lower half-space is intentionally quiet");
verifySubstring(testCase, body, "rigid baffle");
verifySubstring(testCase, body, "high-order impedance boundary is mandatory");
verifySubstring(testCase, body, "do not use or name a Kelvin");
verifySubstring(testCase, body, "reduced FEM ODE");
verifySubstring(testCase, body, "ode45");
verifySubstring(testCase, body, "damping-ratio");
verifySubstring(testCase, body, "not a cavity");
verifySubstring(testCase, body, "decaying membrane/shell vibration");
verifySubstring(testCase, body, "same modeling split can be implemented in NGSolve");
verifySubstring(testCase, body, "retarded boundary");
verifySubstring(testCase, body, "must NOT split the observation field by source direction");
verifySubstring(testCase, body, "all evaluated at every");
verifySubstring(testCase, body, "direction-only painting is a");
verifySubstring(testCase, body, "not a cavity");
verifySubstring(testCase, body, "pressure DOF");
verifySubstring(testCase, body, "lower-half radiation");
verifySubstring(testCase, body, "3D axisymmetric");
verifySubstring(testCase, body, "r-z slice");
verifySubstring(testCase, body, "not yet a full 3D structural-FEM/acoustic-BEM drum mesh");
verifySubstring(testCase, body, "frequency-domain Helmholtz FEM/BEM");
verifySubstring(testCase, body, "inverse FFT");
verifySubstring(testCase, body, "not a periodic sine-wave animation");
verifySubstring(testCase, body, "parallel acoustic-volume teaching lane");
verifySubstring(testCase, body, "not the preferred");
verifySubstring(testCase, body, "BDF generating function");
verifySubstring(testCase, body, "Laplace-domain single-layer");
verifySubstring(testCase, body, "real Lubich CQ TD-BEM");
verifySubstring(testCase, body, "volFemBemCoupledConvolutionQuadrature");
verifySubstring(testCase, body, "H1/P1 interior wave FEM");
verifySubstring(testCase, body, "(1/2 Mb-K(s))*T");
verifySubstring(testCase, body, "-S(s)q + D(s)Tu");
verifySubstring(testCase, body, "Calderon/Johnson-Nedelec coupled CQ");
verifySubstring(testCase, body, "retarded double-layer K(s)");
verifySubstring(testCase, body, "SingleLayerTeaching");
verifySubstring(testCase, body, "GmshPostExport");
verifySubstring(testCase, body, "replace the interior acoustic volume FEM with structural membrane/shell FEM");
verifySubstring(testCase, body, "does not require Gmsh");
end


function testKnowledgeIncludesCurvedVolGeometryTopic(testCase)
body = acoustic_fembem.fembem_knowledge("curved_vol_geometry");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "superparametric");
verifySubstring(testCase, body, "CurvedPanelQuadrature");
verifySubstring(testCase, body, "curvedSingleLayerDirichletSolve");
verifySubstring(testCase, body, "curvedelements");
end


function testKnowledgeIncludesMatlabExecutionPolicy(testCase)
body = acoustic_fembem.fembem_knowledge("matlab_execution_policy");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, ".m functions/scripts");
verifySubstring(testCase, body, "MCP tools");
verifySubstring(testCase, body, "JSON manifests");
end


function testKnowledgeIncludesMathWorksAgenticToolkitPolicy(testCase)
body = acoustic_fembem.fembem_knowledge("mathworks_agentic_toolkit");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "official MathWorks MATLAB MCP Server");
verifySubstring(testCase, body, "MATLAB Agentic Toolkit");
verifySubstring(testCase, body, "runtime");
verifySubstring(testCase, body, "skills");
verifySubstring(testCase, body, "existing-session");
verifySubstring(testCase, body, "acoustic_fembem extension");
verifyEqual(testCase, acoustic_fembem.fembem_knowledge("matlab_mcp_server"), body);
end


function testVolMeshSummaryWrapper(testCase)
out = evalc("acoustic_fembem.check_vol_mesh_summary(""unit_sphere_coarse.vol"")");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_vol_mesh_summary");
verifyEqual(testCase, string(decoded.recommended_gui_viewer), "Netgen/native .vol viewer");
verifyTrue(testCase, contains(string(decoded.recommended_windows_double_click_handler), "Ng_LoadMesh"));
verifyGreaterThan(testCase, decoded.points, 0);
verifyGreaterThan(testCase, decoded.triangles, 0);
verifyGreaterThan(testCase, decoded.tets, 0);
end


function testRepositoryHealthWrapper(testCase)
out = evalc("acoustic_fembem.check_repository_health()");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_repository_health");
verifyEqual(testCase, string(decoded.repository_name), "matlab-acoustic-fembem");
verifyEqual(testCase, decoded.num_validation_cases, 100);
verifyEqual(testCase, decoded.num_verified_cases, 100);
verifyGreaterThanOrEqual(testCase, decoded.num_vol_fixtures, 10);
end


function testResultManifestGateWrapper(testCase)
artifact = completeArtifact();
manifestPath = fullfile(tempdir, "acoustic_fembem_result_manifest_gate_test.json");
fid = fopen(manifestPath, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonencode(artifact));
clear cleanup

out = evalc("acoustic_fembem.check_result_manifest_file(manifestPath, true, true, ""matlab_fem_bem_result_table_v1"")");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_check_result_manifest_file");
end


function testAnalyticAcousticGateWrapper(testCase)
out = evalc("acoustic_fembem.check_fembem_acoustic_gate(""soft"", 2.0, 7, -1)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_acoustic_gate");
verifyEqual(testCase, string(decoded.kind), "soft");
end


function testCrossvalGateWrapper(testCase)
out = evalc("acoustic_fembem.check_fembem_crossval_gate(""galerkin_ngsolve"", ""unit_sphere_coarse.vol"", -1, false)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_crossval_gate");
verifyEqual(testCase, string(decoded.input_format), "netgen_vol_tri_tet");
end


function artifact = completeArtifact()
artifact = struct();
artifact.schema = "matlab-acoustic-fembem.crossval.v1";
artifact.pass = true;
artifact.created_at_utc = "2026-07-04T00:00:00Z";
artifact.versions = struct("matlab", version, "radia_mcp", "test");
artifact.execution = struct( ...
    "run_date_utc", "2026-07-04T00:00:02Z", ...
    "execution_session_id", "MATLAB_TEST");
artifact.expected_created_at_utc = "2026-07-04T00:00:00Z";
artifact.expected_run_date_utc = "2026-07-04T00:00:02Z";
artifact.expected_execution_session_id = "MATLAB_TEST";
artifact.result_output_schema_id = "matlab_fem_bem_result_table_v1";
artifact.result_output_columns = ["alpha", "trace_residual_norm"];
artifact.result_output_units = struct("alpha", "1", "trace_residual_norm", "1");
artifact.timing_breakdown_s = struct("solve", 0.1, "postprocess", 0.02);
artifact.physics_convention_schema_id = "matlab_first_order_fem_bem_coupling_convention_v1";
artifact.postprocess_row_convention_schema_id = "matlab_fem_bem_postprocess_rows_v1";
end
