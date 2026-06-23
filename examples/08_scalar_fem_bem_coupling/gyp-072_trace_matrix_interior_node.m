% GYP-072 trace matrix interior node
% Gypsilab inspiration: doc/FEM-BEM coupling + nonRegressionTest/vibroAcoustic / trace matrix interior node
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-072");
