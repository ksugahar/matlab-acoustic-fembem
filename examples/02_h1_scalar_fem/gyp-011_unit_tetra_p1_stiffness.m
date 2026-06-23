% GYP-011 unit tetra p1 stiffness
% Gypsilab inspiration: openFem + nonRegressionTest/finiteElement / unit tetra p1 stiffness
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-011");
