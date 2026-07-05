repoRoot = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(repoRoot);
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

% Discover every tests/test*.m by globbing the folder, so a newly added test
% cannot silently go unregistered (a manual literal list previously let 7 test
% files - the CQ time-domain, iFFT, drum, and Gmsh-artifact lane - never run).
testEntries = dir(fullfile(repoRoot, "tests", "test*.m"));
testFiles = string(fullfile({testEntries.folder}, {testEntries.name}));

allResults = matlab.unittest.TestResult.empty;
for k = 1:numel(testFiles)
    results = runtests(testFiles(k));
    allResults = [allResults; results(:)]; %#ok<AGROW>
end

disp(table(allResults));
assert(all([allResults.Passed]));
