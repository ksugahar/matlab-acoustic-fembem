% GYP-012 unit tetra p1 mass
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / unit tetra p1 mass
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-012");
