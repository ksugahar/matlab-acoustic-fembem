classdef RwgSpace
%RWGSPACE Lowest-order RWG (triangle-pair) dofs on the boundary surface.
%
%   space = RwgSpace(surface);
%   space.dofEdgeIds                  % interior manifold edges = RWG dofs
%   ids = space.hcurlEdgeIds(hcurl);  % the same dofs inside the volume edge set
%
% Every surface edge is stored (sorted compact node ids); an edge becomes an
% RWG dof only when it is shared by exactly two triangles. The map into the
% volume Nedelec0 edge set is the oriented-edge trace used by H(curl)/RWG
% coupling, so orientation bugs surface before any Maxwell kernel work.

properties (Constant)
    family = "Rao-Wilton-Glisson"
    order = 0
    cell = "triangle-pair"
    basis = "RWG"
end

properties
    surface                % SurfaceMesh carrying compact vtx/tri
    edgesLocal             % all surface edges (nEdges x 2), compact node ids
    edgesGlobal            % the same edges in volume node ids
    triEdges               % edge id per triangle local edge (nTris x 3)
    triEdgeSigns           % +1/-1 local-vs-global orientation (nTris x 3)
    edgeTriangles          % adjacent triangle ids per edge (nEdges x 2, 0 = none)
    oppositeVerticesLocal  % opposite vertex per edge-triangle slot (nEdges x 2)
    dofEdgeIds             % edge ids shared by two triangles (RWG dofs)
    dofEdgesLocal          % dof edges in compact node ids
    dofEdgesGlobal         % dof edges in volume node ids
end

methods
    function space = RwgSpace(surface)
        arguments
            surface (1,1) SurfaceMesh
        end
        space.surface = surface;
        [space.edgesLocal, space.triEdges, space.triEdgeSigns, ...
            space.edgeTriangles, space.oppositeVerticesLocal] = ...
            triEdgeTopology(surface.tri);
        space.edgesGlobal = surface.volNodeIds(space.edgesLocal);

        dofMask = all(space.edgeTriangles > 0, 2);
        space.dofEdgeIds = find(dofMask);
        space.dofEdgesLocal = space.edgesLocal(dofMask, :);
        space.dofEdgesGlobal = space.edgesGlobal(dofMask, :);
    end

    function n = ndof(space)
        n = numel(space.dofEdgeIds);
    end

    function ids = hcurlEdgeIds(space, hcurl)
        %HCURLEDGEIDS Volume Nedelec0 edge id of every RWG dof edge.
        %
        % Errors when a boundary dof edge is missing from the volume edge
        % set; that means the surface does not belong to this volume mesh.
        arguments
            space (1,1) RwgSpace
            hcurl (1,1) Nedelec0Space
        end
        [isTraceEdge, ids] = ismember(space.dofEdgesGlobal, hcurl.edges, "rows");
        if any(~isTraceEdge)
            error("RwgSpace:trace", ...
                "A boundary RWG edge is not present in the volume HCurl edge set.");
        end
    end
end
end


function [edges, triEdges, triEdgeSigns, edgeTriangles, oppositeVertices] = triEdgeTopology(tri)
%TRIEDGETOPOLOGY Unique oriented surface edges with adjacency and signs.

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
        error("RwgSpace:nonmanifold", ...
            "Surface edge belongs to more than two triangles.");
    end
    edgeTriangles(edgeId, slot) = triId;
    oppositeVertices(edgeId, slot) = oppositeRaw(row);
end
end
