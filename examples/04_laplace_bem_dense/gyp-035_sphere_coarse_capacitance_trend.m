% GYP-035 sphere coarse capacitance trend
% Gypsilab inspiration: openOpr + nonRegressionTest/operators + openEbd scalar products / sphere coarse capacitance trend
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-035");
