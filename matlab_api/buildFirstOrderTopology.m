function topology = buildFirstOrderTopology(model)
%BUILDFIRSTORDERTOPOLOGY Build first-order H1/HCurl/RWG connectivity.
%
% Lukas FEM remains a reference. This function creates the minimal clean-room
% topology needed by MATLAB prototypes and radia-ngsolve cross validation:
% volume H1 nodes, volume HCurl edges, boundary scalar P1 nodes, and RWG
% surface-edge dofs.

tet = model.lukas.geo.conn_matrix;
triLocal = model.gypsilab.elt;
traceNodeIds = model.trace.nodeIds(:);

[hcurlEdges, tetEdges, tetEdgeSigns] = buildTetEdges(tet);
[surfaceEdgesLocal, triEdges, triEdgeSigns, edgeTriangles, oppositeVertices] = buildTriEdges(triLocal);

surfaceEdgesGlobal = traceNodeIds(surfaceEdgesLocal);
[isTraceEdge, rwgToHcurlEdgeIds] = ismember(surfaceEdgesGlobal, hcurlEdges, "rows");
if any(~isTraceEdge)
    error("buildFirstOrderTopology:trace", "A boundary RWG edge is not present in the volume HCurl edge set.");
end

rwgDofMask = all(edgeTriangles > 0, 2);

topology = struct();
topology.h1 = struct( ...
    "nodeIds", (1:size(model.lukas.geo.nodes, 1)).', ...
    "traceNodeIds", traceNodeIds);
topology.hcurl = struct( ...
    "edges", hcurlEdges, ...
    "tetEdges", tetEdges, ...
    "tetEdgeSigns", tetEdgeSigns);
topology.scalarBem = struct( ...
    "nodeIds", (1:numel(traceNodeIds)).', ...
    "globalNodeIds", traceNodeIds);
topology.rwg = struct( ...
    "edgesLocal", surfaceEdgesLocal, ...
    "edgesGlobal", surfaceEdgesGlobal, ...
    "triEdges", triEdges, ...
    "triEdgeSigns", triEdgeSigns, ...
    "edgeTriangles", edgeTriangles, ...
    "oppositeVerticesLocal", oppositeVertices, ...
    "dofEdgeIds", find(rwgDofMask), ...
    "dofEdgesLocal", surfaceEdgesLocal(rwgDofMask, :), ...
    "dofEdgesGlobal", surfaceEdgesGlobal(rwgDofMask, :), ...
    "hcurlEdgeIds", rwgToHcurlEdgeIds(rwgDofMask));
topology.trace = struct( ...
    "h1ToScalarBem", sparse(1:numel(traceNodeIds), traceNodeIds, 1, numel(traceNodeIds), size(model.lukas.geo.nodes, 1)), ...
    "rwgToHcurlEdgeIds", rwgToHcurlEdgeIds(rwgDofMask));
topology.policy = "first_order_h1_p1_hcurl_nedelec0_bem_p1_rwg_only";
end


function [edges, tetEdges, tetEdgeSigns] = buildTetEdges(tet)
localPairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
nTets = size(tet, 1);
raw = zeros(6 * nTets, 2);
signsRaw = zeros(6 * nTets, 1);

for e = 1:nTets
    for k = 1:6
        row = (e - 1) * 6 + k;
        pair = tet(e, localPairs(k, :));
        sortedPair = sort(pair);
        raw(row, :) = sortedPair;
        if isequal(pair, sortedPair)
            signsRaw(row) = 1;
        else
            signsRaw(row) = -1;
        end
    end
end

[edges, ~, ic] = unique(raw, "rows");
tetEdges = reshape(ic, 6, nTets).';
tetEdgeSigns = reshape(signsRaw, 6, nTets).';
end


function [edges, triEdges, triEdgeSigns, edgeTriangles, oppositeVertices] = buildTriEdges(tri)
localPairs = [1 2; 2 3; 3 1];
nTri = size(tri, 1);
raw = zeros(3 * nTri, 2);
signsRaw = zeros(3 * nTri, 1);
oppositeRaw = zeros(3 * nTri, 1);

for e = 1:nTri
    for k = 1:3
        row = (e - 1) * 3 + k;
        pair = tri(e, localPairs(k, :));
        sortedPair = sort(pair);
        raw(row, :) = sortedPair;
        if isequal(pair, sortedPair)
            signsRaw(row) = 1;
        else
            signsRaw(row) = -1;
        end
        oppositeRaw(row) = tri(e, setdiff(1:3, localPairs(k, :)));
    end
end

[edges, ~, ic] = unique(raw, "rows");
triEdges = reshape(ic, 3, nTri).';
triEdgeSigns = reshape(signsRaw, 3, nTri).';

edgeTriangles = zeros(size(edges, 1), 2);
oppositeVertices = zeros(size(edges, 1), 2);
for row = 1:numel(ic)
    edgeId = ic(row);
    triId = floor((row - 1) / 3) + 1;
    slot = find(edgeTriangles(edgeId, :) == 0, 1);
    if isempty(slot)
        error("buildFirstOrderTopology:nonmanifold", "Surface edge belongs to more than two triangles.");
    end
    edgeTriangles(edgeId, slot) = triId;
    oppositeVertices(edgeId, slot) = oppositeRaw(row);
end
end
