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


function testKnowledgeIncludesCrossvalTopic(testCase)
body = acoustic_fembem.fembem_knowledge("radia_ngsolve_crossval");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, ".vol");
verifySubstring(testCase, body, "radia-ngsolve");
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
end


function testKnowledgeIncludesMatlabExecutionPolicy(testCase)
body = acoustic_fembem.fembem_knowledge("matlab_execution_policy");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, ".m functions/scripts");
verifySubstring(testCase, body, "MCP tools");
verifySubstring(testCase, body, "JSON manifests");
end


function testVolMeshSummaryWrapper(testCase)
out = evalc("acoustic_fembem.check_vol_mesh_summary(""unit_sphere_coarse.vol"")");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_vol_mesh_summary");
verifyEqual(testCase, string(decoded.recommended_gui_viewer), "Netgen/native .vol viewer");
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
artifact.schema = "cae-ai-lab.crossval.v1";
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
