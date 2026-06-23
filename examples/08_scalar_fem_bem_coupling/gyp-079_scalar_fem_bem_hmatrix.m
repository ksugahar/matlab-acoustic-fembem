% GYP-079 scalar fem bem hmatrix
% Gypsilab inspiration: doc/FEM-BEM coupling + nonRegressionTest/vibroAcoustic / scalar fem bem hmatrix
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-079");
