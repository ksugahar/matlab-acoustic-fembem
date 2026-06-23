% GYP-067 frequency sweep phase
% Gypsilab inspiration: miscellaneous/sphereHelmholtz.m + nonRegressionTest/scattering2d/scattering3d / frequency sweep phase
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-067");
