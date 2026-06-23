% GYP-063 pulsating sphere low ka
% Gypsilab inspiration: miscellaneous/sphereHelmholtz.m + nonRegressionTest/scattering2d/scattering3d / pulsating sphere low ka
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-063");
