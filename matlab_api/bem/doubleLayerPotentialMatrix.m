function rows = doubleLayerPotentialMatrix(surface, points, wavenumber, quadratureOrder)
%DOUBLELAYERPOTENTIALMATRIX D_k[.] at points as a linear map on P1 traces.
%
%   rows = doubleLayerPotentialMatrix(surface, points, k);
%   u_s = rows * t;   % scattered double-layer potential of a P1 trace t
%
% The exterior double-layer potential
%
%   D_k[t](x) = int_S dG_k/dn(y) (x,y) t(y) dS(y)
%
% evaluated at each row of `points`, returned as a (nPoints x nNodes) dense
% matrix so a trace -> field map is available as one linear operator (used by
% rigidScatteringSolve.scatteredAt, the CHIEF interior rows, and the adjoint
% source sensitivity acousticFocusAdjoint). Same split as the operators:
% analytic Laplace panels for the singular part, smooth
% low-frequency-stable correction by quadrature for k > 0, so the k -> 0
% limit is exactly the Laplace evaluation. Outward normals come from the
% mesh orientation signs (fail-loud on unknown orientation).

arguments
    surface (1,1) SurfaceMesh
    points (:,3) double
    wavenumber (1,1) double {mustBeNonnegative}
    quadratureOrder (1,1) double {mustBeMember(quadratureOrder, [1 3 7])} = 3
end

signs = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri;
vtx = surface.vtx;
rows = zeros(size(points, 1), size(vtx, 1));
for t = 1:size(tri, 1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t, :), :), points);
    rows(:, tri(t, :)) = rows(:, tri(t, :)) + signs(t) * J1 / (4 * pi);
end

if wavenumber > 0
    quad = SurfaceQuadrature(surface, quadratureOrder);
    parts = HelmholtzKernel(points, quad.points, ...
        "Wavenumber", wavenumber, ...
        "SourceWeights", quad.weights, ...
        "SourceNormals", quad.outwardNormals());
    rows = rows + parts.doubleLayerSourceNormalCorrection * quad.basis;
end
end
