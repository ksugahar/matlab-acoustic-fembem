% GYP-041 cluster tree split
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / cluster tree split
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-041");
