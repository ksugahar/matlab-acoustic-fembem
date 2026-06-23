% GYP-031 single triangle centroid kernel
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / single triangle centroid kernel
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-031");
