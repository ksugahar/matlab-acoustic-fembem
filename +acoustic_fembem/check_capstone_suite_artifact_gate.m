function out = check_capstone_suite_artifact_gate(summary_json)
%CHECK_CAPSTONE_SUITE_ARTIFACT_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.capstone_suite_artifact_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
