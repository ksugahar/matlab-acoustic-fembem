function out = check_periodic_pm_machine_artifact_gate(summary_json)
%CHECK_PERIODIC_PM_MACHINE_ARTIFACT_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.periodic_pm_machine_artifact_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
