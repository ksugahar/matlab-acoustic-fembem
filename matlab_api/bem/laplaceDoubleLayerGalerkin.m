function K = laplaceDoubleLayerGalerkin(surface, s, soundSpeed, quadratureOrder)
%LAPLACEDOUBLELAYERGALERKIN Laplace-domain double-layer Galerkin BEM matrix K(s).
%
%   K = laplaceDoubleLayerGalerkin(surface, s, soundSpeed, quadratureOrder)
%   K   % (nNodes x nNodes), K_ij = int int phi_i(x) dG/dn_y(x,y;s) phi_j(y)
%
%   The Galerkin double-layer operator for the Laplace-domain (retarded) kernel
%   G(r; s) = exp(-s r/c)/(4 pi r), c = soundSpeed; the sibling of
%   laplaceSingleLayerGalerkin used in the Johnson-Nedelec/Calderon coupling of
%   the CQ solvers (volFemBemCoupled / volFemBemElastic).  The singular Laplace
%   part is integrated analytically per source triangle (laplaceDoubleLayerPanelIntegrals)
%   with the outward-normal orientation sign; the smooth retarded correction is a
%   regular product-Gauss integral whose coincident-point limit is zero.

quad = SurfaceQuadrature(surface, quadratureOrder);
signs = surface.orientation.triangleOrientationSignsToOutward(:);
if any(signs == 0)
    error("laplaceDoubleLayerGalerkin:orientation", ...
        "Surface orientation is unknown for %d triangle(s); cannot assemble K(s).", ...
        sum(signs == 0));
end

nGauss = quad.nPoints();
nNodes = size(surface.vtx, 1);
tri = surface.tri;
vtx = surface.vtx;

P = complex(zeros(nGauss, nNodes));
for t = 1:size(tri, 1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t, :), :), quad.points);
    P(:, tri(t, :)) = P(:, tri(t, :)) + signs(t) * J1;
end
Bw = quad.weightedBasis();
K = Bw.' * P / (4 * pi);

correction = doubleLayerCorrection(quad.points, quad.points, s, ...
    soundSpeed, quad.weights, quad.outwardNormals());
K = K + Bw.' * (correction * quad.basis);
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
