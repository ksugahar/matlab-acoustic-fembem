function h = plotVolMesh(volOrMesh, options)
%PLOTVOLMESH Quick MATLAB figure preview for a first-order .vol mesh.
%
%   plotVolMesh("unit_sphere_coarse.vol")
%   plotVolMesh(mesh, FaceAlpha=0.2)
%
% Netgen remains the preferred interactive viewer for native .vol files.
% This helper is intentionally small: it plots the boundary triangles so a
% MATLAB script or interactive session can confirm scale, orientation, and
% labels before a solve.

arguments
    volOrMesh
    options.FaceAlpha (1,1) double {mustBeGreaterThanOrEqual(options.FaceAlpha, 0), mustBeLessThanOrEqual(options.FaceAlpha, 1)} = 0.35
    options.FaceColor (1,3) double = [0.20 0.55 0.85]
    options.EdgeColor = [0.20 0.20 0.20]
    options.ShowTetEdges (1,1) logical = false
end

mesh = asVolMesh(volOrMesh);

h = patch( ...
    "Faces", mesh.tri, ...
    "Vertices", mesh.vtx, ...
    "FaceColor", options.FaceColor, ...
    "FaceAlpha", options.FaceAlpha, ...
    "EdgeColor", options.EdgeColor);
axis equal
grid on
xlabel("x")
ylabel("y")
zlabel("z")
title(string(mesh.meshId), "Interpreter", "none")
view(3)

if options.ShowTetEdges
    holdState = ishold;
    hold on
    edges = unique(sort([ ...
        mesh.tet(:, [1 2])
        mesh.tet(:, [1 3])
        mesh.tet(:, [1 4])
        mesh.tet(:, [2 3])
        mesh.tet(:, [2 4])
        mesh.tet(:, [3 4])], 2), "rows");
    p = mesh.vtx;
    for k = 1:size(edges, 1)
        pts = p(edges(k, :), :);
        plot3(pts(:, 1), pts(:, 2), pts(:, 3), "-", ...
            "Color", [0.1 0.1 0.1]);
    end
    if ~holdState
        hold off
    end
end
end


function mesh = asVolMesh(value)
if isa(value, "VolMesh")
    mesh = value;
elseif ischar(value) || (isstring(value) && isscalar(value))
    mesh = VolMesh(string(value));
else
    error("plotVolMesh:input", "Input must be a .vol path or a VolMesh object.");
end
end
