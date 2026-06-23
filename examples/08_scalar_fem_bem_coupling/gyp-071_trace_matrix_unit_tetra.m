% GYP-071 trace matrix unit tetra
% Gypsilab inspiration: doc/FEM-BEM coupling + nonRegressionTest/vibroAcoustic / trace matrix unit tetra
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-071");
