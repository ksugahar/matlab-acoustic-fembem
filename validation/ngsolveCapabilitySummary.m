function caps = ngsolveCapabilitySummary()
%NGSOLVECAPABILITYSUMMARY Query NGSolve reference capabilities once.

repoRoot = gypsilabRepoRoot();
script = fullfile(repoRoot, "validation", "ngsolve_capability_summary.py");
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol");
if ~isfolder("C:\temp")
    mkdir("C:\temp");
end
jsonOut = string(tempname("C:\temp")) + ".json";
cleanup = onCleanup(@() deleteIfExists(jsonOut));

cmd = sprintf('python "%s" --vol-file "%s" --output "%s"', script, volFile, jsonOut);
[status, output] = system(cmd);
if status ~= 0
    error("ngsolveCapabilitySummary:python", "NGSolve capability query failed: %s", output);
end
caps = jsondecode(fileread(jsonOut));
if ~caps.ok
    error("ngsolveCapabilitySummary:ngsolve", "NGSolve capability query failed: %s", caps.errorMessage);
end
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
