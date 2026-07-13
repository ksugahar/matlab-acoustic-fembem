function tests = testSlot602GlobalLocalOptimizationGate
tests = functiontests(localfunctions);
end

function testPositive(testCase)
row = baseRow();
result = jsondecode(acoustic_fembem.check_global_local_optimization_replay_gate(jsonencode(row)));
verifyEqual(testCase, string(result.status), "ok");
end

function testNegative(testCase)
row = baseRow(); row.public_gate_status = 'needs_attention'; row.source_copy_preserved = false;
result = jsondecode(acoustic_fembem.check_global_local_optimization_replay_gate(jsonencode(row)));
verifyEqual(testCase, string(result.status), "needs_attention");
verifyFalse(testCase, result.checks.public_gate_passed);
end

function row = baseRow()
row = struct('shared_session','MATLAB_1','matlab_version','release','global_toolbox_available',true, ...
    'optimization_toolbox_available',true,'source_copy_preserved',true, ...
    'source_algorithm_sha256',repmat('a',1,64),'source_objective_sha256',repmat('b',1,64), ...
    'seed_count',3,'short_iterations',220,'long_iterations',600,'public_gate_status','ok', ...
    'run_date_utc','2026-07-11T00:00:00Z');
end
