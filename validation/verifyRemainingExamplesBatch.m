function results = verifyRemainingExamplesBatch(options)
%VERIFYREMAININGEXAMPLESBATCH Verify GYP-011 through GYP-100 in one pass.
%
% This is the second campaign gate after the mesh/topology 10/100 batch. It
% keeps one NGSolve capability snapshot and applies category-specific MATLAB
% teaching checks to the remaining 90 examples.

arguments
    options.WriteLog (1,1) logical = true
end

repoRoot = gypsilabRepoRoot();
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));

cases = validationCatalog();
selected = cases(11:100);
caps = ngsolveCapabilitySummary();
results = repmat(emptyResult(), numel(selected), 1);

for k = 1:numel(selected)
    results(k) = verifyCatalogCase(selected(k), caps);
end

if options.WriteLog
    writeRemainingLog(results, caps);
end

assert(all([results.passed]));
end


function result = emptyResult()
result = struct();
result.id = "";
result.title = "";
result.category = "";
result.reference = "";
result.tolerance = NaN;
result.passed = false;
result.failures = strings(0, 1);
result.details = struct();
end


function writeRemainingLog(results, caps)
logPath = fullfile(tempdir, "gypsilab-validation", ...
    "gypsilab_remaining_90of100_20260624.md");
logDir = fileparts(logPath);
if ~isfolder(logDir)
    mkdir(logDir);
end

lines = [
    "# Gypsilab 100-case campaign: remaining 90/100"
    ""
    "Date: 2026-06-24"
    ""
    "Gate: MATLAB readable FEM/BEM teaching operators against radia-ngsolve/NGSolve capability checks."
    ""
    "Scope:"
    "- GYP-011 through GYP-100"
    "- first-order H1 P1 tetra FEM"
    "- first-order HCurl/Nedelec0 tetra FEM"
    "- scalar P1/RWG trace scaffolds"
    "- dense Laplace and acoustic Helmholtz BEM kernels"
    "- readable Laplace H-matrix blocks and matvecs"
    "- NGSolve.BEM availability smoke gates"
    ""
    "NGSolve snapshot:"
    sprintf("- version: %s", string(caps.version))
    sprintf("- mesh: vertices=%d elements=%d edges=%d", caps.meshVertices, caps.meshElements, caps.meshEdges)
    sprintf("- H1 dofs: %d", caps.h1Dofs)
    sprintf("- HCurl dofs: %d", caps.hcurlDofs)
    sprintf("- BEM: LaplaceSL=%d HelmholtzSL=%d MaxwellSL=%d", caps.hasLaplaceSL, caps.hasHelmholtzSL, caps.hasMaxwellSL)
    ""
    "| ID | Category | Title | Result | Gate |"
    "| --- | --- | --- | --- | --- |"
];

for k = 1:numel(results)
    r = results(k);
    if r.passed
        outcome = "PASS";
    else
        outcome = "FAIL: " + strjoin(r.failures, "; ");
    end
    gate = "";
    if isfield(r.details, "gate")
        gate = string(r.details.gate);
    end
    lines(end + 1, 1) = sprintf("| %s | %s | %s | %s | %s |", ...
        r.id, r.category, r.title, outcome, gate); %#ok<AGROW>
end

lines = [
    lines
    ""
    "Summary:"
    sprintf("- passed: %d / %d", nnz([results.passed]), numel(results))
    "- verified cases: GYP-011 through GYP-100"
    "- GYP-001 through GYP-010 remain covered by `gypsilab_mesh_topology_10of100_20260624.md`."
    "- These are educational verification gates, not performance claims."
    "- The public log records only analytic and open numerical reference gates."
    ""
];

fid = fopen(logPath, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", lines);
clear cleanup
end
