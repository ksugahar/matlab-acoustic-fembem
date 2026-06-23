% GYP-033 laplace symmetry weighted
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / laplace symmetry weighted
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-033");
