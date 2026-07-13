function tests = testCqScatterCausalityMcp
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(repoRoot));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testWrapperRunsCausalCqScatter(testCase)
out = evalc("acoustic_fembem.check_soft_sphere_cq_causality(16, 0.4, 24)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_soft_sphere_cq_causality");
verifyGreaterThanOrEqual(testCase, decoded.result.peak_offset_steps, -1);
verifyLessThanOrEqual(testCase, decoded.result.peak_offset_steps, 3);
verifyLessThan(testCase, decoded.result.max_relative_residual, 1e-6);
end


function testWrapperRejectsUndersampledTimeGrid(testCase)
verifyError(testCase, ...
    @() acoustic_fembem.check_soft_sphere_cq_causality(3, 0.4, 24), ...
    "MATLAB:validators:mustBeGreaterThan");
end
