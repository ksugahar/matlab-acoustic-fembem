function check_cq_response_reality(summary_json, residual_tolerance, imaginary_tolerance)
%CHECK_CQ_RESPONSE_REALITY Print an MCP-facing coupled CQ verdict.
arguments
    summary_json (1,1) string
    residual_tolerance (1,1) double = 1e-10
    imaginary_tolerance (1,1) double = 1e-10
end
report = cqResponseRealityManifest(jsondecode(summary_json), residual_tolerance, imaginary_tolerance);
disp(jsonencode(struct("tool", "acoustic_fembem_cq_response_reality", "ok", report.ok, "result", report)));
end
