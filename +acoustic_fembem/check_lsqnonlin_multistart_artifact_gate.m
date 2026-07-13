function out = check_lsqnonlin_multistart_artifact_gate(summary_json)
%CHECK_LSQNONLIN_MULTISTART_ARTIFACT_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.lsqnonlin_multistart_artifact_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
