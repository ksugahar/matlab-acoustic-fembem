% GYP-060 ngsolve bem low frequency comparison
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / ngsolve bem low frequency comparison
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-060");
