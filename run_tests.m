repoRoot = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

testFiles = [
    fullfile(repoRoot, "tests", "testReadVolTriTet.m")
    fullfile(repoRoot, "tests", "testEducationalMeshImportQuality.m")
    fullfile(repoRoot, "tests", "testFirstOrderFemBemSpaces.m")
    fullfile(repoRoot, "tests", "testEducationalHMatrix.m")
    fullfile(repoRoot, "tests", "testEducationalAcoustics.m")
    fullfile(repoRoot, "tests", "testEducationalCoulombGauge.m")
    fullfile(repoRoot, "tests", "testEducationalOptimization.m")
    fullfile(repoRoot, "tests", "testEducationalGeometricIntegration.m")
    fullfile(repoRoot, "tests", "testValidationCatalog.m")
    fullfile(repoRoot, "tests", "testMeshTopologyValidation.m")
    fullfile(repoRoot, "tests", "testRemainingExamplesValidation.m")
];

allResults = matlab.unittest.TestResult.empty;
for k = 1:numel(testFiles)
    results = runtests(testFiles(k));
    allResults = [allResults; results(:)]; %#ok<AGROW>
end

disp(table(allResults));
assert(all([allResults.Passed]));
