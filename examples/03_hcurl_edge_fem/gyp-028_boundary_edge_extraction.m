% GYP-028 boundary edge extraction
% Gypsilab inspiration: openFem/femNedelec + nonRegressionTest/finiteElement/rtFemRwgNed.m / boundary edge extraction
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-028");
