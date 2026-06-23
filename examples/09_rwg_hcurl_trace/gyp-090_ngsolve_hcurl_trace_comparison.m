% GYP-090 ngsolve hcurl trace comparison
% Gypsilab inspiration: openFem/femRaoWiltonGlisson + femNedelec + femBemDielectrique / ngsolve hcurl trace comparison
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-090");
