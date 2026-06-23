% GYP-099 ngsolve fem bem coupling smoke
% Gypsilab inspiration: Gypsilab nonRegressionTest capstones mirrored against NGSolve.BEM / ngsolve fem bem coupling smoke
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-099");
