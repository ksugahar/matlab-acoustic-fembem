% GYP-058 low frequency sphere monopole
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / low frequency sphere monopole
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-058");
