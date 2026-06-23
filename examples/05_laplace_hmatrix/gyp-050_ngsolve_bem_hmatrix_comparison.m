% GYP-050 ngsolve bem hmatrix comparison
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / ngsolve bem hmatrix comparison
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-050");
