function rows = laplaceDoubleLayerPotential(surface, points, s, soundSpeed, quadratureOrder)
%LAPLACEDOUBLELAYERPOTENTIAL Laplace-domain double-layer potential D(s) at points.
%
%   rows = laplaceDoubleLayerPotential(surface, points, s, soundSpeed, quadratureOrder)
%   rows   % (nPoints x nNodes): rows*p = the outgoing field of a nodal boundary
%   trace p at the (exterior) points, kernel dG/dn_y with G = exp(-s r/c)/(4 pi r).
%
%   The evaluation companion of laplaceDoubleLayerGalerkin: singular Laplace part
%   integrated analytically per source triangle (laplaceDoubleLayerPanelIntegrals)
%   with the outward-normal sign, smooth retarded correction by product Gauss
%   quadrature.  Shared by the CQ solvers (volFemBemCoupled / volFemBemElastic).

signs = surface.orientation.triangleOrientationSignsToOutward(:);
if any(signs == 0)
    error("laplaceDoubleLayerPotential:orientation", ...
        "Surface orientation is unknown for %d triangle(s); cannot evaluate D(s).", ...
        sum(signs == 0));
end

tri = surface.tri;
vtx = surface.vtx;
rows = complex(zeros(size(points, 1), size(vtx, 1)));
for t = 1:size(tri, 1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t, :), :), points);
    rows(:, tri(t, :)) = rows(:, tri(t, :)) + signs(t) * J1 / (4 * pi);
end

quad = SurfaceQuadrature(surface, quadratureOrder);
correction = doubleLayerCorrection(points, quad.points, s, ...
    soundSpeed, quad.weights, quad.outwardNormals());
rows = rows + correction * quad.basis;
end


function C = doubleLayerCorrection(targetPoints, sourcePoints, s, soundSpeed, ...
        sourceWeights, sourceNormals)
%DOUBLELAYERCORRECTION Smooth part of the retarded double-layer kernel, with weights.
nTarget = size(targetPoints, 1);
nSource = size(sourcePoints, 1);
C = complex(zeros(nTarget, nSource));
alpha = s / soundSpeed;
for i = 1:nTarget
    for j = 1:nSource
        delta = targetPoints(i, :) - sourcePoints(j, :);
        r = norm(delta);
        if r == 0
            value = 0.0;      % the double-layer smooth correction vanishes at r=0
        else
            normalDot = dot(delta, sourceNormals(j, :));
            base = normalDot / r^3;
            z = -alpha * r;
            value = base * stableExpTimesOneMinusZMinusOne(z);
        end
        C(i, j) = sourceWeights(j) * value / (4 * pi);
    end
end
end


function value = stableExpTimesOneMinusZMinusOne(z)
%STABLEEXPTIMESONEMINUSZMINUSONE exp(z)*(1-z)-1, Taylor near z=0.
if abs(z) < 1e-5
    value = 0;
    for k = 2:10
        value = value + (1 - k) * z^k / factorial(k);
    end
else
    value = exp(z) * (1 - z) - 1;
end
end
