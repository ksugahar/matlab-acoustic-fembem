% GYP-015 dirichlet elimination
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / dirichlet elimination
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-015");
