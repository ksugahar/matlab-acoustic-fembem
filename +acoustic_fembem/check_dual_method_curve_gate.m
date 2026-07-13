function out = check_dual_method_curve_gate(summary_json)
%CHECK_DUAL_METHOD_CURVE_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.dual_method_curve_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
