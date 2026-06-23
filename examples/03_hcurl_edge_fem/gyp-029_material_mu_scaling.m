% GYP-029 material mu scaling
% Gypsilab inspiration: openFem/femNedelec + nonRegressionTest/finiteElement/rtFemRwgNed.m / material mu scaling
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-029");
