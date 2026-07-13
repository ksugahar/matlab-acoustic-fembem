function out = check_rwg_hcurl_trace_artifact_gate(summary_json)
%CHECK_RWG_HCURL_TRACE_ARTIFACT_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.rwg_hcurl_trace_artifact_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
