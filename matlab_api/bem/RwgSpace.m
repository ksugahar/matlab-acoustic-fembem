classdef RwgSpace
%RWGSPACE Lowest-order RWG (triangle-pair) dofs on the boundary surface.
%
%   space = RwgSpace(surface);
%   space.dofEdgeIds                  % interior manifold edges = RWG dofs
%   ids = space.hcurlEdgeIds(hcurl);  % the same dofs inside the volume edge set
%   [Bx, By, Bz] = space.basisAtQuadrature(quad);
%   [G, detail] = space.gram();       % int f_i . f_j dS (exact, 3-pt rule)
%
% Every surface edge is stored (sorted compact node ids); an edge becomes an
% RWG dof only when it is shared by exactly two triangles. The map into the
% volume Nedelec0 edge set is the oriented-edge trace used by H(curl)/RWG
% coupling, so orientation bugs surface before any Maxwell kernel work.
%
% Basis convention (Rao-Wilton-Glisson 1982): dof edge e with plus triangle
% edgeTriangles(e,1) and minus triangle edgeTriangles(e,2),
%
%   f_e(y) = +l_e/(2 A_plus)  (y - p_plus)    on the plus triangle
%   f_e(y) = -l_e/(2 A_minus) (y - p_minus)   on the minus triangle
%
% where p_+- is the vertex opposite the edge and l_e the edge length, so
% unit current crosses e from plus to minus and div f_e = +-l_e/A_+-.
% The trace identity locked by the tests is
%
%   n x (tangential trace of the volume Nedelec0 edge function)
%     = triEdgeSigns(T,e) * sigma_pm(T,e) * f_e / l_e     on each triangle T
%
% with n the OUTWARD normal, triEdgeSigns the local-winding-vs-ascending
% orientation sign, and sigma_pm = +1/-1 on the plus/minus triangle.

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
    dofEdgeLengths         % length of each dof edge (nDof x 1)
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
        edgeVec = surface.vtx(space.dofEdgesLocal(:, 2), :) ...
            - surface.vtx(space.dofEdgesLocal(:, 1), :);
        space.dofEdgeLengths = sqrt(sum(edgeVec.^2, 2));
    end

    function n = ndof(space)
        n = numel(space.dofEdgeIds);
    end

    function [Bx, By, Bz] = basisAtQuadrature(space, quad)
        %BASISATQUADRATURE RWG basis values at surface quadrature points.
        %
        % Returns three sparse (nPoints x nDof) component matrices so
        % Galerkin operators can sandwich kernels between RWG dofs.
        arguments
            space (1,1) RwgSpace
            quad (1,1) SurfaceQuadrature
        end
        areas = space.surface.areas();
        nPts = quad.nPoints();
        rows = []; cols = []; vx = []; vy = []; vz = [];
        for d = 1:space.ndof()
            e = space.dofEdgeIds(d);
            for slot = 1:2
                t = space.edgeTriangles(e, slot);
                sigma = 3 - 2 * slot;              % +1 on plus, -1 on minus
                pOpp = space.surface.vtx(space.oppositeVerticesLocal(e, slot), :);
                coef = sigma * space.dofEdgeLengths(d) / (2 * areas(t));
                pts = find(quad.triangleIndex == t);
                vals = coef * (quad.points(pts, :) - pOpp);
                rows = [rows; pts(:)]; %#ok<AGROW>
                cols = [cols; repmat(d, numel(pts), 1)]; %#ok<AGROW>
                vx = [vx; vals(:, 1)]; %#ok<AGROW>
                vy = [vy; vals(:, 2)]; %#ok<AGROW>
                vz = [vz; vals(:, 3)]; %#ok<AGROW>
            end
        end
        Bx = sparse(rows, cols, vx, nPts, space.ndof());
        By = sparse(rows, cols, vy, nPts, space.ndof());
        Bz = sparse(rows, cols, vz, nPts, space.ndof());
    end

    function [G, detail] = gram(space)
        %GRAM RWG mass matrix int f_i . f_j dS, exact with the 3-point rule.
        quad = SurfaceQuadrature(space.surface, 3);
        [Bx, By, Bz] = space.basisAtQuadrature(quad);
        W = spdiags(quad.weights, 0, quad.nPoints(), quad.nPoints());
        G = Bx.' * W * Bx + By.' * W * By + Bz.' * W * Bz;
        detail = struct( ...
            "family", space.basis, ...
            "cell", space.cell, ...
            "gram", G, ...
            "quadratureOrder", 3, ...
            "dofEdgeLengths", space.dofEdgeLengths);
    end

    function C = rotatedTraceMap(space, hcurl)
        %ROTATEDTRACEMAP RWG coefficients of the rotated Nedelec0 trace.
        %
        % For a volume field u = sum_E alpha_E N_E the rotated tangential
        % trace n x u|_Gamma (n outward) lies EXACTLY in the RWG space:
        %
        %   n x u|_Gamma = sum_d (C * alpha)_d f_d
        %
        % with C sparse (nDof x hcurl.ndof), C(d, E_d) = gamma_d / l_d and
        % gamma_d = -signOut(T) * triEdgeSigns(T, e) * sigma_pm(T, e)
        % evaluated on either adjacent triangle (tangential continuity
        % makes the two sides agree). This is the H(curl)/RWG coupling
        % contract of the ladder, locked pointwise by the tests.
        arguments
            space (1,1) RwgSpace
            hcurl (1,1) Nedelec0Space
        end
        ids = space.hcurlEdgeIds(hcurl);
        signsOut = space.surface.orientation.triangleOrientationSignsToOutward(:);
        if any(signsOut == 0)
            error("RwgSpace:orientation", ...
                "Surface orientation is unknown; cannot fix the outward normal.");
        end
        gamma = zeros(space.ndof(), 1);
        for d = 1:space.ndof()
            e = space.dofEdgeIds(d);
            t = space.edgeTriangles(e, 1);                 % plus triangle
            le = find(space.triEdges(t, :) == e, 1);
            gamma(d) = -signsOut(t) * space.triEdgeSigns(t, le);   % sigma = +1
        end
        C = sparse(1:space.ndof(), ids, gamma ./ space.dofEdgeLengths, ...
            space.ndof(), hcurl.ndof());
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
