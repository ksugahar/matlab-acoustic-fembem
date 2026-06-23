function summary = ngsolveVolSummary(volFile)
%NGSOLVEVOLSUMMARY Summarize .vol topology through NGSolve.

arguments
    volFile (1,1) string
end

repoRoot = gypsilabRepoRoot();
script = fullfile(repoRoot, "validation", "ngsolve_vol_summary.py");
tempRoot = "C:\temp";
if ~isfolder(tempRoot)
    mkdir(tempRoot);
end
jsonOut = string(tempname(tempRoot)) + ".json";
cleanup = onCleanup(@() deleteIfExists(jsonOut));

cmd = sprintf('python "%s" --output "%s" "%s"', script, jsonOut, volFile);
[status, output] = system(cmd);
if status ~= 0
    error("ngsolveVolSummary:python", "NGSolve summary failed: %s", output);
end

rows = jsondecode(fileread(jsonOut));
summary = rows(1);
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
