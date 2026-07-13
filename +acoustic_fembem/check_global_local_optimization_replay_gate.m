function out = check_global_local_optimization_replay_gate(summary_json)
%CHECK_GLOBAL_LOCAL_OPTIMIZATION_REPLAY_GATE Metadata-only source-tool contract.
try
    row = jsondecode(char(summary_json));
    checks = struct();
    checks.shared_existing_session = startsWith(string(row.shared_session), "MATLAB_");
    checks.release_recorded = strlength(string(row.matlab_version)) > 0;
    checks.global_toolbox_available = logical(row.global_toolbox_available);
    checks.optimization_toolbox_available = logical(row.optimization_toolbox_available);
    checks.source_copy_preserved = logical(row.source_copy_preserved);
    checks.source_algorithm_recorded = strlength(string(row.source_algorithm_sha256)) == 64;
    checks.source_objective_recorded = strlength(string(row.source_objective_sha256)) == 64;
    checks.three_seed_short_long_replay = row.seed_count >= 3 && row.short_iterations < row.long_iterations;
    checks.public_gate_passed = string(row.public_gate_status) == "ok";
    checks.result_date_recorded = strlength(string(row.run_date_utc)) > 0;
    names = fieldnames(checks);
    values = cellfun(@(n) logical(checks.(n)), names);
    issues = names(~values);
    result = struct('schema','acoustic-fembem-global-local-optimization-replay/v1', ...
        'status', ternary(all(values),'ok','needs_attention'), 'checks',checks, ...
        'issues',{issues}, 'notes',{{ ...
        'The legacy multi-swarm search is a basin finder, not proof of a global optimum.', ...
        'Keep the shared MATLAB session, stochastic replay, local polish, and derivative evidence distinct.'}});
catch exception
    result = struct('status','invalid_input','error',exception.message);
end
out = jsonencode(result, PrettyPrint=true);
end

function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end
