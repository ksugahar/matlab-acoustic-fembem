function tests = testCqResponseRealityManifest
tests = functiontests(localfunctions);
end
function testPositiveAndNegative(testCase)
root = fileparts(fileparts(mfilename("fullpath"))); addpath(genpath(fullfile(root,"matlab_api")));
s = struct("method","BDF2","coupling_form","JohnsonNedelec","num_time_steps",16, ...
    "num_laplace_solves",16,"min_real_laplace_parameter",3,"max_relative_residual",4e-18, ...
    "interior_imag_relative",4e-14,"exterior_imag_relative",2e-14,"double_layer_included",true);
good = cqResponseRealityManifest(s); verifyTrue(testCase, good.ok);
s.exterior_imag_relative = 1e-2; bad = cqResponseRealityManifest(s); verifyFalse(testCase, bad.ok);
end
