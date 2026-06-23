% GYP-082 rwg open edge rejection
% Gypsilab inspiration: openFem/femRaoWiltonGlisson + femNedelec + femBemDielectrique / rwg open edge rejection
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-082");
