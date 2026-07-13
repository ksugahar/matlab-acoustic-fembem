function check_cq_time_grid(dt, n_steps, method, contour_samples, radius)
%CHECK_CQ_TIME_GRID Print a compact MCP-facing CQ grid verdict.

arguments
    dt (1,1) double = 0.01
    n_steps (1,1) double = 100
    method (1,1) string = "BDF2"
    contour_samples (1,1) double = 0
    radius (1,1) double = 0
end

report = cqTimeGridManifest(dt, n_steps, method, contour_samples, radius);
payload = struct( ...
    "tool", "acoustic_fembem_cq_time_grid", ...
    "ok", report.ok, ...
    "result", report);
disp(jsonencode(payload));
end
