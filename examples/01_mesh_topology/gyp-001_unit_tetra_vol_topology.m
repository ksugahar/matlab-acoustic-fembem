repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(fullfile(repoRoot, "validation"));
result = runMeshTopologyExample("GYP-001");
