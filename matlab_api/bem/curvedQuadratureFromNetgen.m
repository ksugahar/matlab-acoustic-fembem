function quad = curvedQuadratureFromNetgen(surface, jsonPath, order)
%CURVEDQUADRATUREFROMNETGEN Curved-panel quadrature from a netgen curved-node JSON.
%
%   quad = curvedQuadratureFromNetgen(surface, "ng_sphere_curved_p2.json", 7);
%
% Consumes the convention-free curved-boundary-node companion written by
% tools/export_curved_boundary_nodes.py (NGSolve evaluates each boundary
% triangle's curved geometry map at the Lagrange nodes and stores the physical
% COORDINATES -- no Netgen coefficient basis).  Each exported triangle is
% matched to a SurfaceMesh triangle by its corner coordinates (exact, tolerance
% 1e-7, raises on a miss -- no silent fallback), the curved edge midpoints are
% placed into the [m12 m23 m31] slots by matching corner pairs, and a curve
% order 2 CurvedPanelQuadrature is built on the explicit nodes.  This is the
% optional "netgen .vol -> MATLAB high accuracy" path for GENERAL geometry
% (the self-generated analytic projection only exists for analytic surfaces).

arguments
    surface (1,1) SurfaceMesh
    jsonPath (1,1) string
    order (1,1) double {mustBeMember(order, [1 3 7])} = 7
end

data = jsondecode(fileread(jsonPath));
if data.curve_order ~= 2
    error("curvedQuadratureFromNetgen:order", ...
        "only curve order 2 companions are supported (got %d).", data.curve_order);
end
tris = data.triangles;
nExp = numel(tris);

keys = cell(nExp, 1);
for e = 1:nExp
    keys{e} = cornerKey(tris(e).corners);
end
lookup = containers.Map(keys, 1:nExp);

tol = 1e-7;
nTris = size(surface.tri, 1);
geomNodes = zeros(nTris, 6, 3);
for t = 1:nTris
    C = surface.vtx(surface.tri(t, :), :);       % 3 x 3 MATLAB corner coords
    key = cornerKey(C);
    if ~isKey(lookup, key)
        error("curvedQuadratureFromNetgen:match", ...
            "no exported triangle matches surface triangle %d.", t);
    end
    ex = tris(lookup(key));
    geomNodes(t, 1, :) = C(1, :);
    geomNodes(t, 2, :) = C(2, :);
    geomNodes(t, 3, :) = C(3, :);
    geomNodes(t, 4, :) = edgeMidFor(ex, C(1, :), C(2, :), tol);
    geomNodes(t, 5, :) = edgeMidFor(ex, C(2, :), C(3, :), tol);
    geomNodes(t, 6, :) = edgeMidFor(ex, C(3, :), C(1, :), tol);
end

quad = CurvedPanelQuadrature(surface, order, "GeomNodes", geomNodes);
end


function k = cornerKey(C)
%CORNERKEY Order-independent char key of 3 corner coordinates (sorted rounded rows).
R = sortrows(round(C * 1e7) / 1e7);
k = sprintf('%.7f,', R.');   % char (containers.Map needs char keys)
end


function mid = edgeMidFor(ex, A, B, tol)
%EDGEMIDFOR Curved midpoint of the exported edge joining corners A and B.
% Exported edgemid e joins corners e and mod(e,3)+1 (positional pairing).
for e = 1:3
    a = ex.corners(e, :);
    b = ex.corners(mod(e, 3) + 1, :);
    if (norm(a - A) < tol && norm(b - B) < tol) || ...
       (norm(a - B) < tol && norm(b - A) < tol)
        mid = reshape(ex.edgemids(e, :), 1, 1, 3);
        return
    end
end
error("curvedQuadratureFromNetgen:edge", ...
    "no exported edge matches corner pair (%g %g %g)-(%g %g %g).", A, B);
end
