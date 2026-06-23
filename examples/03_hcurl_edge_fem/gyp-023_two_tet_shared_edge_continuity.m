% GYP-023 two tet shared edge continuity
% Gypsilab inspiration: openFem/femNedelec + nonRegressionTest/finiteElement/rtFemRwgNed.m / two tet shared edge continuity
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-023");
