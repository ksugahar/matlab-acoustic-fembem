% GYP-062 plane wave boundary data
% Gypsilab inspiration: miscellaneous/sphereHelmholtz.m + nonRegressionTest/scattering2d/scattering3d / plane wave boundary data
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-062");
