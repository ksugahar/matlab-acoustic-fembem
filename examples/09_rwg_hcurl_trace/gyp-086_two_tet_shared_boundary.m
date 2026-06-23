% GYP-086 two tet shared boundary
% Gypsilab inspiration: openFem/femRaoWiltonGlisson + femNedelec + femBemDielectrique / two tet shared boundary
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-086");
