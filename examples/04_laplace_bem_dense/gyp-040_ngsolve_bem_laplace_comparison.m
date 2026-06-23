% GYP-040 ngsolve bem laplace comparison
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / ngsolve bem laplace comparison
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-040");
