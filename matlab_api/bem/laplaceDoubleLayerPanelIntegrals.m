function [J0, J1] = laplaceDoubleLayerPanelIntegrals(vertices, points)
%LAPLACEDOUBLELAYERPANELINTEGRALS Analytic double-layer integrals, P1.
%
%   [J0, J1] = laplaceDoubleLayerPanelIntegrals(vertices, points)
%
%   vertices : 3x3, one triangle vertex per row (normal from the CCW order)
%   points   : n x 3 observation points
%   J0       : n x 1,  J0(p)   = int_T  n(y).(x_p - y)/|x_p - y|^3      dS(y)
%   J1       : n x 3,  J1(p,k) = int_T  lambda_k(y) n.(x_p - y)/r^3     dS(y)
%
% For a flat triangle n.(x - y) = w0 (the signed height of x), so the
% kernel is  w0/r^3 = -d(1/r)/dw0  and the closed forms follow by
% differentiating the verified single-layer panel integrals in w0 with the
% in-plane projection held fixed:
%
%   J0 = -dI0/dw0 = signed solid angle of T seen from x
%        (van Oosterom & Strackee 1983, robust arctangent form)
%   d(Ivec)/dw0 = sum_e mHat_e * w0 * f2_e   (the log edge terms)
%   J1_k = lambda_k(rho0) * J0 - w0 * sum_e (grad(lambda_k).mHat_e) f2_e
%
% Partition of unity sum_k J1_k = J0 holds identically. For observation
% points in the triangle plane the principal value is zero and both
% outputs return 0 (the +-1/2 jump lives in the boundary integral
% equation, not in this integral). The kernel carries no 1/(4*pi).

arguments
    vertices (3,3) double
    points (:,3) double
end

p1 = vertices(1, :); p2 = vertices(2, :); p3 = vertices(3, :);
nVec = cross(p2 - p1, p3 - p1);
twoArea = norm(nVec);
if twoArea <= eps
    error("laplaceDoubleLayerPanelIntegrals:degenerate", ...
        "Degenerate source triangle.");
end
nHat = nVec / twoArea;
diam = max([norm(p1 - p2), norm(p2 - p3), norm(p3 - p1)]);

gradLambda = zeros(3, 3);
edgesOpposite = [p3 - p2; p1 - p3; p2 - p1];
for k = 1:3
    gradLambda(k, :) = cross(nHat, edgesOpposite(k, :)) / twoArea;
end

edgeStart = [p1; p2; p3];
edgeEnd = [p2; p3; p1];

nPts = size(points, 1);
J0 = zeros(nPts, 1);
J1 = zeros(nPts, 3);

for p = 1:nPts
    x = points(p, :);
    w0 = dot(x - p1, nHat);
    if abs(w0) <= 1e-12 * diam
        continue   % in-plane principal value is zero
    end
    rho0 = x - w0 * nHat;

    % signed solid angle (van Oosterom & Strackee)
    a = p1 - x; b = p2 - x; c = p3 - x;
    na = norm(a); nb = norm(b); nc = norm(c);
    numerator = det([a; b; c]);
    denominator = na * nb * nc + dot(a, b) * nc + dot(b, c) * na + dot(c, a) * nb;
    J0(p) = -2 * atan2(numerator, denominator);

    edgeTerm = zeros(1, 3);
    for e = 1:3
        aE = edgeStart(e, :);
        bE = edgeEnd(e, :);
        sHat = (bE - aE) / norm(bE - aE);
        mHat = cross(sHat, nHat);
        t0 = dot(aE - rho0, mHat);
        sMinus = dot(aE - rho0, sHat);
        sPlus = dot(bE - rho0, sHat);
        R0sq = t0^2 + w0^2;
        RMinus = sqrt(R0sq + sMinus^2);
        RPlus = sqrt(R0sq + sPlus^2);
        f2 = log((RPlus + sPlus) / (RMinus + sMinus));
        edgeTerm = edgeTerm + mHat * f2;
    end

    for k = 1:3
        lambdaAtRho0 = lambdaAffine(k, rho0, p1, gradLambda);
        J1(p, k) = lambdaAtRho0 * J0(p) ...
            - w0 * dot(gradLambda(k, :), edgeTerm);
    end
end
end


function value = lambdaAffine(k, y, p1, gradLambda)
%LAMBDAAFFINE Affine extension of barycentric lambda_k evaluated at y.

if k == 1
    base = 1.0;
else
    base = 0.0;
end
value = base + dot(gradLambda(k, :), y - p1);
end
