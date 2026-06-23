function tests = testMeshTopologyValidation
%TESTMESHTOPOLOGYVALIDATION Tests for the first verified 10/100 cases.

tests = functiontests(localfunctions);
end


function testMeshTopologyBatchPasses(testCase)
results = verifyMeshTopologyBatch("WriteLog", false);

verifyEqual(testCase, numel(results), 10);
verifyTrue(testCase, all([results.passed]));
verifyEqual(testCase, [results.id].', "GYP-" + compose("%03d", (1:10).'));
end


function testSingleStudentExampleRuns(testCase)
result = runMeshTopologyExample("GYP-001", "Display", false);

verifyTrue(testCase, result.passed);
verifyEqual(testCase, result.matlab.points, result.ngsolve.points);
verifyEqual(testCase, result.matlab.tets, result.ngsolve.tets);
end
