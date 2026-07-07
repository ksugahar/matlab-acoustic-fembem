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
correction = singleLayerCorrection(points, quad.points, s, soundSpeed, quad.weights);
rows = rows + correction * quad.basis;
end


function C = singleLayerCorrection(targetPoints, sourcePoints, s, soundSpeed, sourceWeights)
%SINGLELAYERCORRECTION Smooth part (exp(-s r/c) - 1)/(4 pi r) with weights.
nTarget = size(targetPoints, 1);
nSource = size(sourcePoints, 1);
C = complex(zeros(nTarget, nSource));
alpha = s / soundSpeed;
for i = 1:nTarget
    for j = 1:nSource
        r = norm(targetPoints(i, :) - sourcePoints(j, :));
        if r == 0
            value = -alpha;      % finite limit -alpha/(4 pi) of the smooth part
        else
            z = -alpha * r;
            value = stableExpm1OverR(z, r);
        end
        C(i, j) = sourceWeights(j) * value / (4 * pi);
    end
end
end


function value = stableExpm1OverR(z, r)
%STABLEEXPM1OVERR (exp(z) - 1)/r, Taylor-stable for small |z|.
if abs(z) < 1e-5
    value = 0;
    term = 1;
    for k = 1:10
        term = term * z / k;
        value = value + term / r;
    end
else
    value = (exp(z) - 1) / r;
end
end
