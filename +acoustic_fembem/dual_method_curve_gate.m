function result = dual_method_curve_gate(summary)
%DUAL_METHOD_CURVE_GATE Gate a shared-session two-method curve artifact.

checks = struct();
checks.shared_existing_session = startsWith(string(summary.shared_session), "MATLAB_");
checks.release_recorded = strlength(string(summary.matlab_version)) > 0;
checks.shared_solver_session_survived = logical(summary.shared_solver_available_before) && logical(summary.shared_solver_available_after);
checks.source_kind_recorded = string(summary.source_kind) == "source_native_matlab_example";
sourceDigest = string(summary.source_sha256);
checks.source_digest_recorded = isscalar(sourceDigest) && strlength(sourceDigest) == 64;
checks.source_copy_preserved = logical(summary.source_copy_preserved);
checks.method_ids_recorded = isequal(sort(string(summary.method_ids(:))), ...
    sort(["air_gap_volume_integral"; "stress_contour_integral"]));
checks.sample_count_sufficient = double(summary.sample_count) >= 9;
checks.live_solve_recorded = double(summary.solve_s) > 0;
checks.public_gate_passed = string(summary.public_gate_status) == "ok";
checks.result_date_recorded = strlength(string(summary.run_date_utc)) > 0;

names = fieldnames(checks);
values = cellfun(@(name) logical(checks.(name)), names);
result = struct( ...
    "schema", "matlab-dual-method-curve-artifact/v1", ...
    "status", ternary(all(values), "ok", "needs_attention"), ...
    "checks", checks, ...
    "issues", {names(~values)}, ...
    "notes", {{ ...
        "Use the existing shared MATLAB session; do not start a second process for a source validation.", ...
        "Two-method agreement is reusable for torque, force, acoustic power, and dense-versus-compressed operator curves.", ...
        "Preserve method identifiers, units, source digest, and the independent public gate result."}});
end


function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end
