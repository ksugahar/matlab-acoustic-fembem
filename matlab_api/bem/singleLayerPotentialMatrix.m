function rows = singleLayerPotentialMatrix(surface, points, wavenumber, quadratureOrder)
%SINGLELAYERPOTENTIALMATRIX S_k[.] at points as a linear map on P1 densities.
%
%   rows = singleLayerPotentialMatrix(surface, points, k);
%   u = rows * q;      % single-layer potential of a P1 boundary density q
%
% The exterior single-layer potential
%
%   S_k[q](x) = int_S G_k(x,y) q(y) dS(y),   G_k = exp(1i*k*r)/(4*pi*r)
%
% evaluated at each row of `points`, returned as a (nPoints x nNodes) dense
% matrix (the density -> field map as one linear operator). Same split as
% the operators: the singular Laplace kernel is integrated ANALYTICALLY over
% every source triangle (laplacePanelIntegrals) and the smooth
% low-frequency-stable correction by quadrature for k > 0, so the k -> 0
% limit is exactly the Laplace evaluation. Companion of
% doubleLayerPotentialMatrix; used by singleLayerDirichletSolve and the FSI
% exterior representation p_s = D_k[p] - S_k[q].

arguments
    surface (1,1) SurfaceMesh
    points (:,3) double
    wavenumber (1,1) double {mustBeNonnegative}
    quadratureOrder (1,1) double {mustBeMember(quadratureOrder, [1 3 7])} = 3
end

tri = surface.tri;
vtx = surface.vtx;
rows = zeros(size(points, 1), size(vtx, 1));
for t = 1:size(tri, 1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), points);
    rows(:, tri(t, :)) = rows(:, tri(t, :)) + I1 / (4 * pi);
end

if wavenumber > 0
    quad = SurfaceQuadrature(surface, quadratureOrder);
    parts = HelmholtzKernel(points, quad.points, ...
        "Wavenumber", wavenumber, "SourceWeights", quad.weights);
    rows = rows + parts.singleLayerCorrection * quad.basis;
end
end
