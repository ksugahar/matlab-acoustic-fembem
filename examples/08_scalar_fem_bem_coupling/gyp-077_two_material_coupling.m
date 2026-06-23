% GYP-077 two material coupling
% Gypsilab inspiration: doc/FEM-BEM coupling + nonRegressionTest/vibroAcoustic / two material coupling
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-077");
