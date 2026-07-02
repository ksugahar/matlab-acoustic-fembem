classdef SurfaceMesh
%SURFACEMESH Compact boundary triangle mesh of a VolMesh.
%
% Boundary P1 dofs use compact node ids 1..nNodes, and volNodeIds maps each
% compact node back to its one-based volume node id:
%
%   surface = mesh.boundary();       % or SurfaceMesh(mesh)
%   surface.vtx(k, :)                % coordinates of compact boundary node k
%   surface.volNodeIds(k)            % the volume node it came from
%   surface.areas()                  % triangle areas
%
% Each triangle row keeps its boundary-condition identity (boundary number,
% boundary name, adjacent tetrahedron) so BEM kernels and coupling notebooks
% can name the surface rows they use.

properties
    vtx           % compact boundary node coordinates (nNodes x 3)
    tri           % triangles in compact node ids (nTris x 3)
    triGlobal     % the same triangles in volume node ids (nTris x 3)
    col           % Netgen boundary number per triangle (nTris x 1)
    volNodeIds    % compact node id -> one-based volume node id (nNodes x 1)
    names         % boundary name per triangle (nTris x 1 string)
    rowIdentity   % per-triangle identity: nodes, boundary number/name, adjacent tet
    orientation   % boundary orientation summary inherited from the VolMesh
    surfaceMeshId % mesh identity of this boundary trace mesh
end

methods
    function surface = SurfaceMesh(mesh)
        arguments
            mesh (1,1) VolMesh
        end
        volNodeIds = unique(mesh.tri(:));
        globalToLocal = zeros(max(volNodeIds), 1);
        globalToLocal(volNodeIds) = 1:numel(volNodeIds);

        surface.vtx = mesh.vtx(volNodeIds, :);
        surface.tri = globalToLocal(mesh.tri);
        surface.triGlobal = mesh.tri;
        surface.col = mesh.triCol;
        surface.volNodeIds = volNodeIds;
        surface.names = mesh.boundaryNames();
        surface.rowIdentity = boundaryRowIdentity( ...
            mesh.tri, ...
            mesh.triCol, ...
            surface.names, ...
            mesh.boundaryOrientation.adjacentTetIndices);
        surface.orientation = mesh.boundaryOrientation;
        surface.surfaceMeshId = mesh.meshId + ":boundary_tri_p1";
    end

    function a = areas(surface)
        %AREAS Triangle areas 0.5 * norm((v2 - v1) x (v3 - v1)).
        v1 = surface.vtx(surface.tri(:, 1), :);
        v2 = surface.vtx(surface.tri(:, 2), :);
        v3 = surface.vtx(surface.tri(:, 3), :);
        n = cross(v2 - v1, v3 - v1, 2);
        a = 0.5 * sqrt(sum(n.^2, 2));
    end

    function c = centroids(surface)
        %CENTROIDS Triangle centroids, the collocation points of this lane.
        c = (surface.vtx(surface.tri(:, 1), :) ...
            + surface.vtx(surface.tri(:, 2), :) ...
            + surface.vtx(surface.tri(:, 3), :)) / 3;
    end

    function mesh = gypsilabMsh(surface)
        %GYPSILABMSH Real Gypsilab msh(vtx, elt, col) for interop studies.
        %
        % Requires the Gypsilab openMsh folder on the MATLAB path; errors
        % loudly instead of returning a placeholder.
        if exist("msh", "file") ~= 2
            error("SurfaceMesh:gypsilab", ...
                "Gypsilab msh class is not on the path; addpath <gypsilab>/openMsh first.");
        end
        mesh = msh(surface.vtx, surface.tri, surface.col);
    end
end
end


function identity = boundaryRowIdentity(triGlobal, boundaryNumbers, boundaryNames, adjacentTetIndices)
%BOUNDARYROWIDENTITY Bind each boundary triangle row to its BC labels.

nRows = size(triGlobal, 1);
identity = repmat(struct( ...
    "surface_triangle_index", 0, ...
    "surface_triangle_nodes", zeros(1, 3), ...
    "boundary_number", 0, ...
    "boundary_name", "", ...
    "adjacent_tet_index", 0), nRows, 1);
for k = 1:nRows
    identity(k).surface_triangle_index = k;
    identity(k).surface_triangle_nodes = triGlobal(k, :);
    identity(k).boundary_number = boundaryNumbers(k);
    identity(k).boundary_name = string(boundaryNames(k));
    identity(k).adjacent_tet_index = adjacentTetIndices(k);
end
end
