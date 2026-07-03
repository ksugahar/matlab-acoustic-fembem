function [I0, I1] = laplacePanelIntegrals(vertices, points)
%LAPLACEPANELINTEGRALS Analytic 1/r and P1/r integrals over one triangle.
%
%   [I0, I1] = laplacePanelIntegrals(vertices, points)
%
%   vertices : 3x3, one triangle vertex per row
%   points   : n x 3 observation points
%   I0       : n x 1,  I0(p)   = int_T  1 / |x_p - y|      dS(y)
%   I1       : n x 3,  I1(p,k) = int_T  lambda_k(y) / |x_p - y| dS(y)
%
% Closed forms after Wilton et al. (1984): per-edge log terms plus the
% signed-height arctangent term for the constant integral, and the in-plane
% edge-normal vector integral for the linear part.  With the barycentric
% function written as an affine field,
%
%   lambda_k(y) = lambda_k(rho0) + grad(lambda_k) . (y - rho0),
%
% the linear integral is  I1(:,k) = lambda_k(rho0) * I0 + grad(lambda_k) . Ivec
% where Ivec = int_T (y - rho0)/|x - y| dS and rho0 is the projection of the
% observation point onto the triangle plane.
%
% The kernel is 1/r (no 1/(4*pi)); the caller owns the Green constant.
% Observation points may sit anywhere, including inside the triangle
% (the 1/r surface integral is finite there); points exactly ON an edge
% line of the triangle are the only excluded case and raise an error.

arguments
    vertices (3,3) double
    points (:,3) double
end

p1 = vertices(1, :); p2 = vertices(2, :); p3 = vertices(3, :);
nVec = cross(p2 - p1, p3 - p1);
twoArea = norm(nVec);
if twoArea <= eps
    error("laplacePanelIntegrals:degenerate", "Degenerate source triangle.");
end
nHat = nVec / twoArea;

% barycentric gradients in the triangle plane (constant, one per vertex)
gradLambda = zeros(3, 3);
edgesOpposite = [p3 - p2; p1 - p3; p2 - p1];   % edge opposite to vertex k
for k = 1:3
    gradLambda(k, :) = cross(nHat, edgesOpposite(k, :)) / twoArea;
end

% counter-clockwise edge loop (with respect to nHat)
edgeStart = [p1; p2; p3];
edgeEnd = [p2; p3; p1];

nPts = size(points, 1);
I0 = zeros(nPts, 1);
I1 = zeros(nPts, 3);

for p = 1:nPts
    x = points(p, :);
    w0 = dot(x - p1, nHat);        % signed height above the plane
    rho0 = x - w0 * nHat;          % in-plane projection

    logSum = 0.0;
    betaSum = 0.0;
    vecSum = zeros(1, 3);

    for e = 1:3
        a = edgeStart(e, :);
        b = edgeEnd(e, :);
        sHat = (b - a) / norm(b - a);
        mHat = cross(sHat, nHat);              % outward in-plane edge normal
        t0 = dot(a - rho0, mHat);              % signed distance to edge line
        sMinus = dot(a - rho0, sHat);
        sPlus = dot(b - rho0, sHat);
        R0sq = t0^2 + w0^2;
        RMinus = sqrt(R0sq + sMinus^2);        % = |x - a|
        RPlus = sqrt(R0sq + sPlus^2);          % = |x - b|

        if RMinus + RPlus <= norm(b - a) * (1 + 1e-12)
            error("laplacePanelIntegrals:onEdge", ...
                "Observation point lies on an edge line of the triangle.");
        end
        % log((RPlus + sPlus)/(RMinus + sMinus)), stabilized when s < 0
        f2 = log((RPlus + sPlus) / (RMinus + sMinus));

        logSum = logSum + t0 * f2;
        betaSum = betaSum ...
            + atan2(t0 * sPlus, R0sq + abs(w0) * RPlus) ...
            - atan2(t0 * sMinus, R0sq + abs(w0) * RMinus);
        vecSum = vecSum + mHat * 0.5 * (R0sq * f2 + sPlus * RPlus - sMinus * RMinus);
    end

    I0(p) = logSum - abs(w0) * betaSum;
    for k = 1:3
        lambdaAtRho0 = lambdaAffine(k, rho0, p1, gradLambda);
        I1(p, k) = lambdaAtRho0 * I0(p) + dot(gradLambda(k, :), vecSum);
    end
end
end


function value = lambdaAffine(k, y, p1, gradLambda)
%LAMBDAAFFINE Affine extension of barycentric lambda_k evaluated at y.

% lambda_k(p1) is 1 for k == 1 and 0 otherwise; extend affinely from p1.
if k == 1
    base = 1.0;
else
    base = 0.0;
end
value = base + dot(gradLambda(k, :), y - p1);
end
