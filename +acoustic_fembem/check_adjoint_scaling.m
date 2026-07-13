function check_adjoint_scaling(rows_json, gradient_tolerance, affine_tolerance, minimum_final_ratio)
%CHECK_ADJOINT_SCALING Print an MCP-facing reverse-mode scaling verdict.
arguments
    rows_json (1,1) string
    gradient_tolerance (1,1) double = 1e-6
    affine_tolerance (1,1) double = 1e-10
    minimum_final_ratio (1,1) double = 1.0
end
rows = jsondecode(rows_json);
report = adjointScalingManifest(rows, gradient_tolerance, affine_tolerance, minimum_final_ratio);
disp(jsonencode(struct("tool", "acoustic_fembem_adjoint_scaling", "ok", report.ok, "result", report)));
end
