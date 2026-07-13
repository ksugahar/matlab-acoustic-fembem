function results = verifyMeshTopologyBatch(options)
%VERIFYMESHTOPOLOGYBATCH Verify the first 10 mesh/topology examples.

arguments
    options.WriteLog (1,1) logical = true
end

repoRoot = gypsilabRepoRoot();
addpath(genpath(fullfile(repoRoot, "matlab_api")));
cases = meshTopologyCaseTable();
results = repmat(emptyResult(), numel(cases), 1);

for k = 1:numel(cases)
    results(k) = verifySingleMeshTopology(cases(k));
end

if options.WriteLog
    writeMeshTopologyLog(results);
end

assert(all([results.passed]));
end


function result = emptyResult()
result = struct();
result.id = "";
result.title = "";
result.volFile = "";
result.expectOk = false;
result.matlab = struct();
result.ngsolve = struct();
result.failures = strings(0, 1);
result.passed = false;
end


function writeMeshTopologyLog(results)
logPath = fullfile(tempdir, "gypsilab-validation", ...
    "gypsilab_mesh_topology_10of100_20260624.md");
logDir = fileparts(logPath);
if ~isfolder(logDir)
    mkdir(logDir);
end
lines = [
    "# Gypsilab 100-case campaign: mesh topology 10/100"
    ""
    "Date: 2026-06-24"
    ""
    "Gate: MATLAB `readVolTriTet` / `FemBemModel` vs NGSolve `Mesh(.vol)`."
    ""
    "| ID | Title | Expected | MATLAB | NGSolve | Result |"
    "| --- | --- | --- | --- | --- | --- |"
];
for k = 1:numel(results)
    r = results(k);
    lines(end + 1, 1) = sprintf("| %s | %s | %s | %s | %s | %s |", ...
        r.id, r.title, okText(r.expectOk), okText(r.matlab.ok), okText(r.ngsolve.ok), passText(r.passed)); %#ok<AGROW>
end
lines = [
    lines
    ""
    "Summary:"
    sprintf("- passed: %d / %d", nnz([results.passed]), numel(results))
    "- verified cases: GYP-001 through GYP-010"
    "- negative cases GYP-008/GYP-009 are verified as MATLAB tri/tet policy rejections."
    "- NGSolve may accept broader Netgen element families; the education solver intentionally does not."
    ""
];
fid = fopen(logPath, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", lines);
clear cleanup
end


function text = okText(tf)
if tf
    text = "ok";
else
    text = "reject";
end
end


function text = passText(tf)
if tf
    text = "PASS";
else
    text = "FAIL";
end
end
