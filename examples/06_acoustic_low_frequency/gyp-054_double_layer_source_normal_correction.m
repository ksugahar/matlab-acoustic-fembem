% GYP-054 double layer source normal correction
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / double layer source normal correction
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-054");
