% GYP-038 normal derivative sign
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / normal derivative sign
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-038");
