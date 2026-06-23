% GYP-053 single layer taylor correction
% Gypsilab inspiration: openEbd Helmholtz kernels + radiationImpedances + acoustic papers / single layer taylor correction
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-053");
