function tests = testMcpDualMethodCurveTool
%TESTMCPDUALMETHODCURVETOOL MCP extension and artifact gate tests.
tests = functiontests(localfunctions);
end


function testExtensionManifestNamesJsonWrapper(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
path = fullfile(repoRoot, "mcp", "extensions", "matlab-dual-method-curve-tools.json");
manifest = jsondecode(fileread(path));
verifyEqual(testCase, string(manifest.tools.name), "matlab_dual_method_curve_artifact_gate");
verifyEqual(testCase, string(manifest.signatures.matlab_dual_method_curve_artifact_gate.function), ...
    "acoustic_fembem.check_dual_method_curve_gate");
end


function testGateAcceptsSharedSessionCurve(testCase)
result = acoustic_fembem.dual_method_curve_gate(validSummary());
verifyEqual(testCase, string(result.status), "ok");
verifyTrue(testCase, result.checks.shared_solver_session_survived);
end


function testJsonWrapperRejectsLostSessionAndMissingMethod(testCase)
summary = validSummary();
summary.shared_solver_available_after = false;
summary.method_ids = "air_gap_volume_integral";
result = jsondecode(acoustic_fembem.check_dual_method_curve_gate(jsonencode(summary)));
verifyEqual(testCase, string(result.status), "needs_attention");
verifyFalse(testCase, result.checks.shared_solver_session_survived);
verifyFalse(testCase, result.checks.method_ids_recorded);
end


function summary = validSummary()
summary = struct( ...
    "shared_session", "MATLAB_1212", ...
    "matlab_version", "R2026a", ...
    "shared_solver_available_before", true, ...
    "shared_solver_available_after", true, ...
    "source_kind", "source_native_matlab_example", ...
    "source_sha256", repmat('a', 1, 64), ...
    "source_copy_preserved", true, ...
    "method_ids", ["air_gap_volume_integral", "stress_contour_integral"], ...
    "sample_count", 21, ...
    "solve_s", 35, ...
    "public_gate_status", "ok", ...
    "run_date_utc", "2026-07-12T00:00:00Z");
end
