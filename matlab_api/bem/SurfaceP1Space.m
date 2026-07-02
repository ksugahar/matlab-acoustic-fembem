classdef SurfaceP1Space
%SURFACEP1SPACE First-order scalar P1 space on the boundary triangles.
%
%   space = SurfaceP1Space(surface);
%   [M, detail] = space.mass();   % int u v dS with the exact P1 rule
%
% The dofs are the compact boundary nodes of the SurfaceMesh. The singular
% Laplace kernel is deliberately assembled elsewhere (HMatrix and the
% acoustic operators); this class owns the regular surface mass only.

properties (Constant)
    family = "Lagrange"
    order = 1
    cell = "triangle"
    basis = "P1"
end

properties
    surface   % SurfaceMesh carrying compact vtx/tri
end

methods
    function space = SurfaceP1Space(surface)
        arguments
            surface (1,1) SurfaceMesh
        end
        space.surface = surface;
    end

    function n = ndof(space)
        n = size(space.surface.vtx, 1);
    end

    function [M, detail] = mass(space)
        %MASS Assemble int u v dS with the exact P1 triangle mass rule.
        vtx = space.surface.vtx;
        tri = space.surface.tri;
        nNodes = size(vtx, 1);
        nTri = size(tri, 1);

        ii = zeros(9 * nTri, 1);
        jj = zeros(9 * nTri, 1);
        vv = zeros(9 * nTri, 1);
        areas = zeros(nTri, 1);
        normals = zeros(nTri, 3);
        local = zeros(3, 3, nTri);

        cursor = 1;
        for e = 1:nTri
            ids = tri(e, :);
            X = vtx(ids, :);
            [Me, area, normal] = localP1TriMass(X);
            local(:, :, e) = Me;
            areas(e) = area;
            normals(e, :) = normal;

            [I, J] = ndgrid(ids, ids);
            span = cursor:(cursor + 8);
            ii(span) = I(:);
            jj(span) = J(:);
            vv(span) = Me(:);
            cursor = cursor + 9;
        end

        M = sparse(ii, jj, vv, nNodes, nNodes);
        detail = struct( ...
            "family", space.basis, ...
            "cell", space.cell, ...
            "surfaceMass", M, ...
            "localMass", local, ...
            "areas", areas, ...
            "normals", normals, ...
            "totalArea", sum(areas), ...
            "globalNodeIds", space.surface.volNodeIds);
    end
end
end


function [Me, area, normal] = localP1TriMass(X)
%LOCALP1TRIMASS Exact P1 mass matrix for one triangle.

crossVec = cross(X(2, :) - X(1, :), X(3, :) - X(1, :));
normCross = norm(crossVec);
area = 0.5 * normCross;
if area <= eps
    error("SurfaceP1Space:degenerate", "Degenerate triangle with near-zero area.");
end

normal = crossVec ./ normCross;
Me = area / 12 * [2 1 1; 1 2 1; 1 1 2];
end
