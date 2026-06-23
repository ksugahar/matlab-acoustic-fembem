repoRoot = fileparts(mfilename("fullpath"));
addpath(fullfile(repoRoot, "matlab_api"));

testFiles = [
    fullfile(repoRoot, "tests", "testReadVolTriTet.m")
    fullfile(repoRoot, "tests", "testFirstOrderFemBemSpaces.m")
    fullfile(repoRoot, "tests", "testEducationalHMatrix.m")
    fullfile(repoRoot, "tests", "testEducationalAcoustics.m")
];

allResults = matlab.unittest.TestResult.empty;
for k = 1:numel(testFiles)
    results = runtests(testFiles(k));
    allResults = [allResults; results(:)]; %#ok<AGROW>
end

disp(table(allResults));
assert(all([allResults.Passed]));
