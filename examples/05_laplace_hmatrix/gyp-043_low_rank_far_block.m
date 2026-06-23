% GYP-043 low rank far block
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / low rank far block
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-043");
