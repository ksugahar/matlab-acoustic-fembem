function check_hmatrix_scaling(point_counts_json, rank_tolerance, matvec_tolerance, max_rank, max_storage_exponent)
%CHECK_HMATRIX_SCALING Run and print an MCP-facing ACA+ scaling verdict.
arguments
    point_counts_json (1,1) string = "[60,120,240]"
    rank_tolerance (1,1) double = 1e-8
    matvec_tolerance (1,1) double = 1e-8
    max_rank (1,1) double = 20
    max_storage_exponent (1,1) double = 1.25
end
pointCounts = reshape(double(jsondecode(point_counts_json)), 1, []);
study = hmatrixScalingStudy(pointCounts, "RankTolerance", rank_tolerance);
report = hmatrixScalingManifest(study.rows, matvec_tolerance, max_rank, max_storage_exponent);
disp(jsonencode(struct( ...
    "tool", "acoustic_fembem_hmatrix_scaling", ...
    "ok", report.ok, ...
    "study", study, ...
    "result", report)));
end
