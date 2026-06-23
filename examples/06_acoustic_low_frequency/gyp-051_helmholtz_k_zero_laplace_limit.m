% GYP-051 helmholtz k zero laplace limit
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / helmholtz k zero laplace limit
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-051");
