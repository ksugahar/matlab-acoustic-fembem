repoRoot = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

testFiles = [
    fullfile(repoRoot, "tests", "testReadVolTriTet.m")
    fullfile(repoRoot, "tests", "testMeshImportQuality.m")
    fullfile(repoRoot, "tests", "testFirstOrderFemBemSpaces.m")
    fullfile(repoRoot, "tests", "testLaplaceDirichletSolve.m")
    fullfile(repoRoot, "tests", "testLaplacePanelIntegrals.m")
    fullfile(repoRoot, "tests", "testGalerkinSingleLayer.m")
    fullfile(repoRoot, "tests", "testFemBemCoupledSolve.m")
    fullfile(repoRoot, "tests", "testFemBemHelmholtzCoupling.m")
    fullfile(repoRoot, "tests", "testNgsolveBemCrossCheck.m")
    fullfile(repoRoot, "tests", "testHelmholtzScattering.m")
    fullfile(repoRoot, "tests", "testSonicCrystalChain.m")
    fullfile(repoRoot, "tests", "testDuctBandGap.m")
    fullfile(repoRoot, "tests", "testRwgVectorCoupling.m")
    fullfile(repoRoot, "tests", "testHMatrix.m")
    fullfile(repoRoot, "tests", "testAcoustics.m")
    fullfile(repoRoot, "tests", "testCoulombGauge.m")
    fullfile(repoRoot, "tests", "testOptimizationGates.m")
    fullfile(repoRoot, "tests", "testGeometricIntegration.m")
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
