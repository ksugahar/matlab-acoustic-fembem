% GYP-064 pulsating sphere mid ka
% Gypsilab inspiration: miscellaneous/sphereHelmholtz.m + nonRegressionTest/scattering2d/scattering3d / pulsating sphere mid ka
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-064");
