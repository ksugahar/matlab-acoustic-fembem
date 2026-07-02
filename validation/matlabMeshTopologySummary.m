function summary = matlabMeshTopologySummary(volFile)
%MATLABMESHTOPOLOGYSUMMARY Summarize .vol topology through the MATLAB API.

arguments
    volFile (1,1) string
end

summary = struct();
summary.source = volFile;
try
    mesh = readVolTriTet(volFile);
    model = FemBemModel(volFile);
    summary.ok = true;
    summary.errorIdentifier = "";
    summary.errorMessage = "";
    summary.points = mesh.summary.points;
    summary.triangles = mesh.summary.triangles;
    summary.tets = mesh.summary.tets;
    summary.materials = mesh.summary.materials;
    summary.boundaries = mesh.summary.boundaries;
    summary.traceNodeCount = numel(mesh.traceNodeIds);
    summary.hcurlEdges = size(model.hcurl.edges, 1);
    summary.rwgDofs = numel(model.rwg.dofEdgeIds);
    summary.materialNames = mapValues(mesh.materials);
    summary.boundaryNames = mapValues(mesh.boundaries);
catch err
    summary.ok = false;
    summary.errorIdentifier = string(err.identifier);
    summary.errorMessage = string(err.message);
    summary.points = NaN;
    summary.triangles = NaN;
    summary.tets = NaN;
    summary.materials = NaN;
    summary.boundaries = NaN;
    summary.traceNodeCount = NaN;
    summary.hcurlEdges = NaN;
    summary.rwgDofs = NaN;
    summary.materialNames = strings(0, 1);
    summary.boundaryNames = strings(0, 1);
end
end


function vals = mapValues(map)
keysCell = keys(map);
if isempty(keysCell)
    vals = strings(0, 1);
    return
end
keysNumeric = cell2mat(keysCell);
keysNumeric = sort(keysNumeric(:));
vals = strings(numel(keysNumeric), 1);
for k = 1:numel(keysNumeric)
    vals(k) = string(map(keysNumeric(k)));
end
end
