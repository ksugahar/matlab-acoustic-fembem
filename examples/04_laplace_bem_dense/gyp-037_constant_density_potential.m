% GYP-037 constant density potential
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / constant density potential
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-037");
