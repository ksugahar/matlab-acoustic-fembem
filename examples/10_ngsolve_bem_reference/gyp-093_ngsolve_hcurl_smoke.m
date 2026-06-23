% GYP-093 ngsolve hcurl smoke
% Gypsilab inspiration: Gypsilab nonRegressionTest capstones mirrored against NGSolve.BEM / ngsolve hcurl smoke
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-093");
