% GYP-048 storage ratio report
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / storage ratio report
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-048");
