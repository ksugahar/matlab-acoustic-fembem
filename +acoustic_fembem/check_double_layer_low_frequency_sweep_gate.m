function out = check_double_layer_low_frequency_sweep_gate(summary_json)
%CHECK_DOUBLE_LAYER_LOW_FREQUENCY_SWEEP_GATE JSON MCP wrapper.
try
    summary = jsondecode(char(summary_json));
    result = acoustic_fembem.double_layer_low_frequency_sweep_gate(summary);
catch exception
    result = struct("status", "invalid_input", "error", exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end
