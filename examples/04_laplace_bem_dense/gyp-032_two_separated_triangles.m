% GYP-032 two separated triangles
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / two separated triangles
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-032");
