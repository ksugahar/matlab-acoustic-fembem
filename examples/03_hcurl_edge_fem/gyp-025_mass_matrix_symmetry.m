% GYP-025 mass matrix symmetry
% Gypsilab inspiration: openFem/femNedelec + nonRegressionTest/finiteElement/rtFemRwgNed.m / mass matrix symmetry
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-025");
