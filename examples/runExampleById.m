function result = runExampleById(caseId, options)
%RUNEXAMPLEBYID Run a verified example or show a planned example card.

arguments
    caseId (1,1) string
    options.Display (1,1) logical = true
end

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "matlab_api"));
addpath(fullfile(repoRoot, "validation"));

cases = validationCatalog();
idx = find([cases.id] == caseId, 1);
if isempty(idx)
    error("runExampleById:case", "Unknown example id: %s", caseId);
end
item = cases(idx);

if item.status == "verified" && startsWith(item.category, "01_mesh_topology")
    result = runMeshTopologyExample(caseId, "Display", options.Display);
    return
end

if item.status == "verified"
    result = verifyCatalogCase(item, ngsolveCapabilitySummary());
    if options.Display
        fprintf("%s %s\n", result.id, result.title);
        fprintf("  status: verified\n");
        fprintf("  category: %s\n", result.category);
        fprintf("  gate: %s\n", result.details.gate);
        fprintf("  passed: %d\n", result.passed);
    end
    return
end

result = struct();
result.id = item.id;
result.title = item.title;
result.status = item.status;
result.category = item.category;
result.gypsilabInspiration = item.gypsilabInspiration;
result.reference = item.reference;
result.secondaryReference = item.secondaryReference;
result.tolerance = item.tolerance;
result.passed = item.status == "verified";

if options.Display
    fprintf("%s %s\n", result.id, result.title);
    fprintf("  status: %s\n", result.status);
    fprintf("  Gypsilab inspiration: %s\n", result.gypsilabInspiration);
    fprintf("  reference gate: %s\n", result.reference);
    if result.secondaryReference ~= ""
        fprintf("  secondary reference: %s\n", result.secondaryReference);
    end
end
end
