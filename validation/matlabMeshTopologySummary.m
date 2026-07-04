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
    summary.materials = countLabels(mesh.tetMat);
    summary.boundaries = countLabels(mesh.triCol);
    summary.traceNodeCount = numel(mesh.traceNodeIds);
    summary.hcurlEdges = size(model.hcurl.edges, 1);
    summary.rwgDofs = numel(model.rwg.dofEdgeIds);
    summary.materialNames = labelNames(mesh.materials, mesh.tetMat);
    summary.boundaryNames = labelNames(mesh.boundaries, mesh.triCol);
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


function n = countLabels(ids)
ids = unique(ids(:));
ids = ids(ids > 0);
n = numel(ids);
end


function vals = labelNames(map, ids)
ids = unique(ids(:));
ids = ids(ids > 0);
if isempty(ids)
    vals = strings(0, 1);
    return
end
vals = strings(numel(ids), 1);
for k = 1:numel(ids)
    id = ids(k);
    if isKey(map, id)
        vals(k) = string(map(id));
    else
        vals(k) = "default";
    end
end
end
