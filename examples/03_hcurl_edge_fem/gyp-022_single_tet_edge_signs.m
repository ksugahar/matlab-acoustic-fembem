% GYP-022 single tet edge signs
% Gypsilab inspiration: openFem/femNedelec + nonRegressionTest/finiteElement/rtFemRwgNed.m / single tet edge signs
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-022");
