% GYP-045 hmatrix matvec vs dense
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / hmatrix matvec vs dense
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-045");
