% GYP-047 leaf size sweep
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / leaf size sweep
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-047");
