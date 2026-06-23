% GYP-056 complex symmetry weighted
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / complex symmetry weighted
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-056");
