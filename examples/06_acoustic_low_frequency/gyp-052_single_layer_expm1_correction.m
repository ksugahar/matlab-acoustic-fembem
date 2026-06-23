% GYP-052 single layer expm1 correction
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / single layer expm1 correction
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-052");
