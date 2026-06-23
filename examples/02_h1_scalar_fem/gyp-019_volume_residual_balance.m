% GYP-019 volume residual balance
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / volume residual balance
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-019");
