% GYP-039 dense matvec reproducibility
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / dense matvec reproducibility
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-039");
