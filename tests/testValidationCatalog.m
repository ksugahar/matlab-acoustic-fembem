function tests = testValidationCatalog
%TESTVALIDATIONCATALOG Tests for the 100-case validation catalog.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));
end


function testCatalogHasOneHundredUniqueCases(testCase)
cases = validationCatalog();

verifyEqual(testCase, numel(cases), 100);
verifyEqual(testCase, numel(unique([cases.id])), 100);
verifyTrue(testCase, all(startsWith([cases.id], "GYP-")));
end


function testEachCategoryHasTenCases(testCase)
cases = validationCatalog();
categories = unique([cases.category]);

verifyEqual(testCase, numel(categories), 10);
for k = 1:numel(categories)
    verifyEqual(testCase, nnz([cases.category] == categories(k)), 10, ...
        "Category " + categories(k) + " should have 10 cases.");
end
end


function testNoUnverifiedCasePretendsToHaveLog(testCase)
cases = validationCatalog();
statuses = [cases.status];

verifyTrue(testCase, all(ismember(statuses, ["planned", "verified", "retired"])));
for k = 1:numel(cases)
    if cases(k).status ~= "verified"
        verifyEqual(testCase, cases(k).validationLog, "");
    end
end
end


function testVerifiedCasesHaveExampleAndLog(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
cases = validationCatalog();
verified = cases([cases.status] == "verified");

verifyEqual(testCase, numel(verified), 100);
for k = 1:numel(verified)
    verifyTrue(testCase, isfile(fullfile(repoRoot, verified(k).examplePath)), ...
        "Missing example for " + verified(k).id);
    verifyNotEqual(testCase, verified(k).validationLog, "");
end


function testAllCasesHaveExampleScriptAndGypsilabInspiration(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
cases = validationCatalog();

for k = 1:numel(cases)
    verifyTrue(testCase, isfile(fullfile(repoRoot, cases(k).examplePath)), ...
        "Missing example for " + cases(k).id);
    verifyNotEqual(testCase, cases(k).gypsilabInspiration, "", ...
        "Missing Gypsilab inspiration for " + cases(k).id);
end
end
end


function testCatalogUsesRadiaNgsolveReference(testCase)
cases = validationCatalog();

verifyTrue(testCase, all([cases.reference] == "radia-ngsolve"));
verifyTrue(testCase, all([cases.tolerance] > 0));
verifyTrue(testCase, all(endsWith([cases.examplePath], ".m")));
end


function testAcousticCasesCarryPublicSecondaryReference(testCase)
cases = validationCatalog();
categories = [cases.category];
secondary = [cases.secondaryReference];
isAcousticReference = ismember(categories, [
    "06_acoustic_low_frequency"
    "07_acoustic_helmholtz"
    "10_ngsolve_bem_reference"
]);

verifyTrue(testCase, all(contains(secondary(isAcousticReference), "open numerical acoustic")));
verifyTrue(testCase, all(secondary(~isAcousticReference) == ""));
end
