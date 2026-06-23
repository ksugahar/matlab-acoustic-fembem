% GYP-094 ngsolve laplace bem smoke
% Gypsilab inspiration: Gypsilab nonRegressionTest capstones mirrored against NGSolve.BEM / ngsolve laplace bem smoke
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-094");
