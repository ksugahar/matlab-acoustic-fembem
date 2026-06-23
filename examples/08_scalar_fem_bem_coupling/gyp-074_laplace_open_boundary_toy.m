% GYP-074 laplace open boundary toy
% Gypsilab inspiration: doc/FEM-BEM coupling + nonRegressionTest/vibroAcoustic / laplace open boundary toy
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-074");
