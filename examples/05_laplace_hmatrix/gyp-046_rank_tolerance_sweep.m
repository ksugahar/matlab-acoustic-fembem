% GYP-046 rank tolerance sweep
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / rank tolerance sweep
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-046");
