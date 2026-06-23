% GYP-059 low frequency hmatrix candidate
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / low frequency hmatrix candidate
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-059");
