% GYP-055 small kr sweep
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / small kr sweep
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-055");
