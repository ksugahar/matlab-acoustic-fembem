% GYP-020 ngsolve p1 scalar comparison
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / ngsolve p1 scalar comparison
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-020");
