% GYP-049 sphere hmatrix laplace
% Gypsilab inspiration: openHmx + nonRegressionTest/hierarchicalMatrix / sphere hmatrix laplace
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-049");
