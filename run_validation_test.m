repoRoot = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(repoRoot);
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

% Heavy validation lane.  Policy (adopted from Radia): tests/ holds the fast
% implementation-regression checks that run every dev loop, and validation_test/
% holds the heavy numerical validation -- physics gates (resonance, radiation
% force, spectral cross-checks, sonic crystal), FEM/BEM coupling, and the CQ
% drum-roll / scattering movies.  Run this before a release; run run_tests.m for
% the quick loop.  Both glob their folder so a new test cannot go unregistered.
testEntries = dir(fullfile(repoRoot, "validation_test", "test*.m"));
testFiles = string(fullfile({testEntries.folder}, {testEntries.name}));

allResults = matlab.unittest.TestResult.empty;
for k = 1:numel(testFiles)
    results = runtests(testFiles(k));
    allResults = [allResults; results(:)]; %#ok<AGROW>
end

disp(table(allResults));
assert(all([allResults.Passed]));
