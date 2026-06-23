% GYP-013 constant gradient patch
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / constant gradient patch
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-013");
