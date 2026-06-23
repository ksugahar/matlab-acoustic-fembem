% GYP-044 dense near block
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / dense near block
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-044");
