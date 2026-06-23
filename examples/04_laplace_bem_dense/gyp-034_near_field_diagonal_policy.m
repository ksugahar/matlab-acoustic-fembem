% GYP-034 near field diagonal policy
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / near field diagonal policy
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-034");
