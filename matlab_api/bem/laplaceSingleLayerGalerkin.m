function V = laplaceSingleLayerGalerkin(surface, s, soundSpeed, quadratureOrder)
%LAPLACESINGLELAYERGALERKIN Laplace-domain single-layer Galerkin BEM matrix V(s).
%
%   V = laplaceSingleLayerGalerkin(surface, s, soundSpeed, quadratureOrder)
%   V   % (nNodes x nNodes), V_ij = int int phi_i(x) G(x,y;s) phi_j(y)
%
% The Galerkin single-layer operator for the Laplace-domain (retarded) kernel
%
%   G(r; s) = exp(-s r / c) / (4 pi r),     c = soundSpeed,
%
% the building block of the Lubich convolution-quadrature time-domain BEM
% (volTdBemConvolutionQuadrature evaluates this at the CQ nodes s = delta(zeta)/dt).
%
% It is the SAME operator as the frequency-domain Helmholtz single layer
% GalerkinSingleLayer: on the imaginary axis s = -1i c k the kernel
% exp(-s r/c) = exp(1i k r), so
%
%   laplaceSingleLayerGalerkin(surface, -1i c k, c, q)
%       == GalerkinSingleLayer(surface, "Wavenumber", k, "QuadratureOrder", q).matrix
%
% to machine precision (locked by tests/testLaplaceSingleLayerGalerkin).  This
% is the CQ correctness anchor: the retarded kernel scaling s/c, the 1/(4 pi),
% and the exponent sign are all pinned to the analytically-validated
% frequency-domain single layer, not just to a self-consistent residual.
%
% Built exactly like GalerkinSingleLayer: the singular Laplace part 1/(4 pi r)
% is integrated ANALYTICALLY over every source triangle (laplacePanelIntegrals)
% at each test Gauss point; the smooth correction (exp(-s r/c) - 1)/(4 pi r) is
% regular and taken by product Gauss quadrature, with the coincident-point value
% left at 0 to match the validated frequency-domain convention (HelmholtzKernel).

arguments
    surface (1,1) SurfaceMesh
    s (1,1) double
    soundSpeed (1,1) double {mustBePositive} = 1.0
    quadratureOrder (1,1) double {mustBeMember(quadratureOrder, [1 3 7])} = 1
end

quad = SurfaceQuadrature(surface, quadratureOrder);
nNodes = size(surface.vtx, 1);
tri = surface.tri;
vtx = surface.vtx;

% --- singular Laplace part: analytic source integral per triangle ---
P = zeros(quad.nPoints(), nNodes);
for t = 1:size(tri, 1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), quad.points);
    P(:, tri(t, :)) = P(:, tri(t, :)) + I1;
end
Bw = quad.weightedBasis();
V = Bw.' * P / (4 * pi);

% --- smooth retarded correction: regular product Gauss quadrature ---
correction = laplaceSingleLayerCorrection(quad.points, quad.points, s, ...
    soundSpeed, quad.weights);
V = V + Bw.' * (correction * quad.basis);
end


function C = laplaceSingleLayerCorrection(targetPoints, sourcePoints, s, soundSpeed, sourceWeights)
%LAPLACESINGLELAYERCORRECTION Smooth part (exp(-s r/c) - 1)/(4 pi r) with weights.
nTarget = size(targetPoints, 1);
nSource = size(sourcePoints, 1);
C = complex(zeros(nTarget, nSource));
alpha = s / soundSpeed;
for i = 1:nTarget
    for j = 1:nSource
        r = norm(targetPoints(i, :) - sourcePoints(j, :));
        if r == 0
            % Coincident quadrature point: leave the smooth correction at 0,
            % the same convention the frequency-domain HelmholtzKernel uses, so
            % V(s) is identical to GalerkinSingleLayer on the imaginary axis.
            value = 0;
        else
            value = stableExpm1OverR(-alpha * r, r);
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
