% GYP-057 low frequency matvec
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / low frequency matvec
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-057");
