function tests = testRemainingExamplesValidation
%TESTREMAININGEXAMPLESVALIDATION Tests for GYP-011 through GYP-100.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));
end


function testRemainingBatchPasses(testCase)
results = verifyRemainingExamplesBatch("WriteLog", false);

verifyEqual(testCase, numel(results), 90);
verifyTrue(testCase, all([results.passed]));
verifyEqual(testCase, [results.id].', "GYP-" + compose("%03d", (11:100).'));
end


function testVerifiedNonMeshExampleRuns(testCase)
result = runExampleById("GYP-011", "Display", false);

verifyTrue(testCase, result.passed);
verifyEqual(testCase, result.category, "02_h1_scalar_fem");
end
