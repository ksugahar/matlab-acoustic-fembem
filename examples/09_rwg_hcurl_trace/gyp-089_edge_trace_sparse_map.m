% GYP-089 edge trace sparse map
% Gypsilab inspiration: openFem/femRaoWiltonGlisson + femNedelec + femBemDielectrique / edge trace sparse map
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-089");
