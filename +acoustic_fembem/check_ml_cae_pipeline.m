function check_ml_cae_pipeline(stage, seed)
%CHECK_ML_CAE_PIPELINE MCP-facing JSON wrapper for the CAE/ML pipeline.
report = acoustic_fembem.ml_cae_pipeline(string(stage), seed);
disp(jsonencode(struct("tool", "acoustic_fembem_ml_cae_pipeline", ...
    "ok", report.ok, "result", report)));
if ~report.ok
    error("acoustic_fembem:MlCaeGateFailed", ...
        "ML/CAE stage '%s' did not pass its gate.", string(stage));
end
end
