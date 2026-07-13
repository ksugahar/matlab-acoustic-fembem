function tests = testMcpDoubleLayerLowFrequencyTool
%TESTMCPDOUBLELAYERLOWFREQUENCYTOOL MCP extension and metadata gate tests.
tests = functiontests(localfunctions);
end


function testExtensionManifestNamesJsonWrapper(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
path = fullfile(repoRoot, "mcp", "extensions", ...
    "acoustic-fembem-double-layer-low-frequency-tools.json");
manifest = jsondecode(fileread(path));
verifyEqual(testCase, string(manifest.tools.name), ...
    "acoustic_fembem_double_layer_low_frequency_sweep_gate");
verifyEqual(testCase, string(manifest.signatures.acoustic_fembem_double_layer_low_frequency_sweep_gate.function), ...
    "acoustic_fembem.check_double_layer_low_frequency_sweep_gate");
end


function testGateAcceptsSharedSessionSweep(testCase)
summary = validSummary();
result = acoustic_fembem.double_layer_low_frequency_sweep_gate(summary);
verifyEqual(testCase, string(result.status), "ok");
verifyTrue(testCase, result.checks.source_examples_pass);
verifyTrue(testCase, result.checks.zero_to_small_kr_covered);
end


function testJsonWrapperRejectsMissingSourceExample(testCase)
summary = validSummary();
summary.source_example_ids = ["GYP-052", "GYP-053"];
result = jsondecode(acoustic_fembem.check_double_layer_low_frequency_sweep_gate(jsonencode(summary)));
verifyEqual(testCase, string(result.status), "needs_attention");
verifyFalse(testCase, result.checks.source_examples_identified);
end


function summary = validSummary()
summary = struct( ...
    "shared_session", "MATLAB_1212", ...
    "matlab_version", "R2026a", ...
    "source_copy_preserved", true, ...
    "source_example_ids", ["GYP-052", "GYP-053", "GYP-054"], ...
    "source_example_pass", [true, true, true], ...
    "kernel_family", "helmholtz_source_normal_double_layer", ...
    "time_convention", "exp(+i*k*r)", ...
    "row_count", 9, ...
    "minimum_kr_abs", 0, ...
    "maximum_kr_abs", 7.5e-4, ...
    "public_gate_status", "ok", ...
    "run_date_utc", "2026-07-12T00:00:00Z");
end
