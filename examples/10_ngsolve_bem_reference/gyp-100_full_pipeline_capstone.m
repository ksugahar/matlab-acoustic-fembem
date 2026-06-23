% GYP-100 full pipeline capstone
% Gypsilab inspiration: Gypsilab nonRegressionTest capstones mirrored against NGSolve.BEM / full pipeline capstone
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-100");
