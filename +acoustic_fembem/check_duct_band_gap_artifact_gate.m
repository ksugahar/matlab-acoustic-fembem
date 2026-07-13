function json = check_duct_band_gap_artifact_gate(summary_json)
%CHECK_DUCT_BAND_GAP_ARTIFACT_GATE JSON entry point for MCP.
arguments
    summary_json (1,1) string
end
summary = jsondecode(summary_json);
json = jsonencode(acoustic_fembem.duct_band_gap_artifact_gate(summary), PrettyPrint=true);
end
