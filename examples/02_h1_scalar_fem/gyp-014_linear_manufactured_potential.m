% GYP-014 linear manufactured potential
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / linear manufactured potential
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-014");
