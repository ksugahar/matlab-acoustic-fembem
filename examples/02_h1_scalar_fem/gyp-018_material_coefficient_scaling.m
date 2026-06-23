% GYP-018 material coefficient scaling
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / material coefficient scaling
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-018");
