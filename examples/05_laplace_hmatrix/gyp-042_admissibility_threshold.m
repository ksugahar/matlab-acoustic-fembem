% GYP-042 admissibility threshold
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / admissibility threshold
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-042");
