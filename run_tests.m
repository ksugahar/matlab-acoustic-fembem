repoRoot = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(repoRoot);
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

% Fast lane.  Policy (adopted from Radia): tests/ holds the fast
% implementation-regression checks that run every dev loop; the heavy numerical
% validation (physics gates, FEM/BEM coupling, drum-roll / scattering movies)
% lives in validation_test/ -- see run_validation_test.m.  Glob every
% tests/test*.m so a newly added test cannot silently go unregistered.
testEntries = dir(fullfile(repoRoot, "tests", "test*.m"));
testFiles = string(fullfile({testEntries.folder}, {testEntries.name}));

allResults = matlab.unittest.TestResult.empty;
for k = 1:numel(testFiles)
    results = runtests(testFiles(k));
    allResults = [allResults; results(:)]; %#ok<AGROW>
end

disp(table(allResults));
assert(all([allResults.Passed]));
