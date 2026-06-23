% GYP-017 two tet patch continuity
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / two tet patch continuity
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-017");
