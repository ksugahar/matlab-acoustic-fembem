function result = simplex_stationarity_artifact_gate(summary)
%SIMPLEX_STATIONARITY_ARTIFACT_GATE Gate a shared-session optimizer audit.

checks = struct();
checks.shared_existing_session = startsWith(string(summary.shared_session), "MATLAB_");
checks.release_recorded = strlength(string(summary.matlab_version)) > 0;
checks.shared_solver_session_survived = logical(summary.shared_solver_available_before) && logical(summary.shared_solver_available_after);
checks.source_kind_recorded = string(summary.source_kind) == "source_native_matlab_optimization_example";
digests = string(struct2cell(summary.source_digests));
checks.source_digests_recorded = numel(digests) >= 3 && all(strlength(digests) == 64);
checks.source_copy_preserved = logical(summary.source_copy_preserved);
checks.live_execution_recorded = double(summary.solve_s) > 0;
checks.finite_difference_step_recorded = double(summary.finite_difference_step) > 0;
checks.method_roles_recorded = isequal(sort(string(summary.method_ids(:))), ...
    sort(["pathological_initial_simplex"; "independent_control"]));
checks.public_gate_passed = string(summary.public_gate_status) == "ok";
checks.false_convergence_detected = logical(summary.public_false_convergence_detected);
checks.independent_control_accepted = logical(summary.public_independent_control_accepted);
checks.result_date_recorded = strlength(string(summary.run_date_utc)) > 0;

names = fieldnames(checks);
values = cellfun(@(name) logical(checks.(name)), names);
result = struct( ...
    "schema", "matlab-simplex-stationarity-artifact/v1", ...
    "status", ternary(all(values), "ok", "needs_attention"), ...
    "checks", checks, ...
    "issues", {names(~values)}, ...
    "notes", {{ ...
        "A small simplex or objective spread is not a stationarity certificate.", ...
        "Use an independently evaluated gradient and at least one accepted control route.", ...
        "Keep the shared solver-connected MATLAB session alive and preserve source digests."}});
end


function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end
