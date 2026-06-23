function bem = assembleGypsilabP1SurfaceMass(model)
%ASSEMBLEGYPSILABP1SURFACEMASS Assemble P1 triangle surface mass.
%
% This is the first-order BEM surface-space check. The singular Laplace kernel
% is deliberately assembled later through Gypsilab regularize/integral.

vtx = model.gypsilab.vtx;
tri = model.gypsilab.elt;
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

bem = struct();
bem.family = "P1";
bem.cell = "triangle";
bem.surfaceMass = sparse(ii, jj, vv, nNodes, nNodes);
bem.localMass = local;
bem.areas = areas;
bem.normals = normals;
bem.totalArea = sum(areas);
bem.globalNodeIds = model.gypsilab.globalNodeIds;
end


function [Me, area, normal] = localP1TriMass(X)
%LOCALP1TRIMASS Exact P1 mass matrix for one triangle.

crossVec = cross(X(2, :) - X(1, :), X(3, :) - X(1, :));
normCross = norm(crossVec);
area = 0.5 * normCross;
if area <= eps
    error("assembleGypsilabP1SurfaceMass:degenerate", "Degenerate triangle with near-zero area.");
end

normal = crossVec ./ normCross;
Me = area / 12 * [2 1 1; 1 2 1; 1 1 2];
end
