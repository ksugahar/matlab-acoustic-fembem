% GYP-016 neumann face load
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / neumann face load
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-016");
