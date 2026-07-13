function result = double_layer_low_frequency_sweep_gate(summary)
%DOUBLE_LAYER_LOW_FREQUENCY_SWEEP_GATE Gate shared-session source evidence.

requiredIds = ["GYP-052", "GYP-053", "GYP-054"];
checks = struct();
checks.shared_existing_session = startsWith(string(summary.shared_session), "MATLAB_");
checks.release_recorded = strlength(string(summary.matlab_version)) > 0;
checks.source_copy_preserved = logical(summary.source_copy_preserved);
observedIds = sort(string(summary.source_example_ids(:)));
checks.source_examples_identified = isequal(observedIds, sort(requiredIds(:)));
checks.source_examples_pass = all(logical(summary.source_example_pass));
checks.kernel_family_recorded = string(summary.kernel_family) == "helmholtz_source_normal_double_layer";
checks.time_convention_recorded = string(summary.time_convention) == "exp(+i*k*r)";
checks.nine_row_sweep_recorded = double(summary.row_count) == 9;
checks.zero_to_small_kr_covered = double(summary.minimum_kr_abs) == 0 && double(summary.maximum_kr_abs) <= 1e-3;
checks.public_gate_passed = string(summary.public_gate_status) == "ok";
checks.result_date_recorded = strlength(string(summary.run_date_utc)) > 0;

names = fieldnames(checks);
values = cellfun(@(name) logical(checks.(name)), names);
result = struct( ...
    "schema", "acoustic-fembem-double-layer-low-frequency-sweep/v1", ...
    "status", ternary(all(values), "ok", "needs_attention"), ...
    "checks", checks, ...
    "issues", {names(~values)}, ...
    "notes", {{ ...
        "The source-normal double-layer regular correction starts at O((k*r)^2).", ...
        "A total-kernel match alone is insufficient because the Laplace term can hide correction cancellation.", ...
        "The source-normal orientation must travel with every reusable result."}});
end


function value = ternary(condition, yesValue, noValue)
if condition, value = yesValue; else, value = noValue; end
end
