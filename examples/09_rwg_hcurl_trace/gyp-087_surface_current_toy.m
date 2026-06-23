% GYP-087 surface current toy
% Gypsilab inspiration: openFem/femRaoWiltonGlisson + femNedelec + femBemDielectrique / surface current toy
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "examples"));
result = runExampleById("GYP-087");
