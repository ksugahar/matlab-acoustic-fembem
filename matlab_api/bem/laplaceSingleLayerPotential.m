function rows = laplaceSingleLayerPotential(surface, points, s, soundSpeed, quadratureOrder)
%LAPLACESINGLELAYERPOTENTIAL Laplace-domain single-layer potential S(s) at points.
%
%   rows = laplaceSingleLayerPotential(surface, points, s, soundSpeed, quadratureOrder)
%   rows   % (nPoints x nNodes): rows*q = the outgoing field of a nodal density q
%   at the (exterior) points, kernel exp(-s r/c)/(4 pi r), c = soundSpeed.
%
%   The evaluation companion of laplaceSingleLayerGalerkin: the singular Laplace
%   part 1/(4 pi r) is integrated ANALYTICALLY over every source triangle
%   (laplacePanelIntegrals) at each target point, and the smooth retarded
%   correction (exp(-s r/c) - 1)/(4 pi r) is taken by product Gauss quadrature.
%   Shared by the CQ solvers (volTdBem / volFemBemCoupled / volFemBemElastic).

tri = surface.tri; vtx = surface.vtx;
rows = zeros(size(points, 1), size(vtx, 1));
for t = 1:size(tri, 1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), points);
    rows(:, tri(t, :)) = rows(:, tri(t, :)) + I1 / (4 * pi);
end

quad = SurfaceQuadrature(surface, quadratureOrder);
correction = singleLayerSmoothCorrection(points, quad.points, s, soundSpeed, quad.weights);
rows = rows + correction * quad.basis;
end
