function out = check_energy_budget_trace_artifact_gate(summary_json)
%CHECK_ENERGY_BUDGET_TRACE_ARTIFACT_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.energy_budget_trace_artifact_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
