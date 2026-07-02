function result = runMeshTopologyExample(caseId, options)
%RUNMESHTOPOLOGYEXAMPLE Run one mesh/topology validation example.

arguments
    caseId (1,1) string
    options.Display (1,1) logical = true
end

repoRoot = gypsilabRepoRoot();
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "validation"));

cases = meshTopologyCaseTable();
idx = find([cases.id] == caseId, 1);
if isempty(idx)
    error("runMeshTopologyExample:case", "Unknown mesh topology case: %s", caseId);
end

result = verifySingleMeshTopology(cases(idx));
if options.Display
    fprintf("%s %s: %s\n", result.id, result.title, passText(result.passed));
    fprintf("  MATLAB ok:  %d\n", result.matlab.ok);
    fprintf("  NGSolve ok: %d\n", result.ngsolve.ok);
end
assert(result.passed);
end


function text = passText(tf)
if tf
    text = "PASS";
else
    text = "FAIL";
end
end
