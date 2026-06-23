% GYP-070 ngsolve bem helmholtz comparison
% Gypsilab inspiration: miscellaneous/sphereHelmholtz.m + nonRegressionTest/scattering2d/scattering3d / ngsolve bem helmholtz comparison
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-070");
