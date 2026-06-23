% GYP-095 ngsolve helmholtz bem smoke
% Gypsilab inspiration: Gypsilab nonRegressionTest capstones mirrored against NGSolve.BEM / ngsolve helmholtz bem smoke
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-095");
