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
% exp(-s r/c) = exp(1i k r).  The two agree exactly EXCEPT at coincident
% quadrature points: this operator keeps the finite limit of the smooth
% correction there (the correct product-Gauss sample, see below), while
% GalerkinSingleLayer / HelmholtzKernel drop it.  The difference is the known
% term
%
%   Delta_ij = (-alpha / (4 pi)) * sum_g w_g^2 phi_i(x_g) phi_j(x_g),   alpha = s/c,
%
% so, locked by tests/testLaplaceSingleLayerGalerkin,
%
%   laplaceSingleLayerGalerkin(surface, -1i c k, c, q)
%       == GalerkinSingleLayer(surface, "Wavenumber", k, "QuadratureOrder", q).matrix + Delta
%
% to machine precision.  This pins the retarded kernel scaling s/c, the 1/(4 pi),
% and the exponent sign to the analytically-validated frequency-domain single
% layer, not just to a self-consistent CQ residual.
%
% Built exactly like GalerkinSingleLayer: the singular Laplace part 1/(4 pi r)
% is integrated ANALYTICALLY over every source triangle (laplacePanelIntegrals)
% at each test Gauss point; the smooth correction (exp(-s r/c) - 1)/(4 pi r) is a
% regular integrand taken by product Gauss quadrature, INCLUDING its finite limit
% -alpha/(4 pi) at coincident points (the correct sample, and the same
% convention as the coupled CQ lane volFemBemCoupledConvolutionQuadrature).

arguments
    surface (1,1) SurfaceMesh
    s (1,1) double {mustBeFinite}
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
            % Coincident quadrature point: the smooth correction
            % (exp(-alpha r) - 1)/(4 pi r) is regular with the finite limit
            % -alpha/(4 pi) as r -> 0.  Keep that limit (the correct product-
            % Gauss sample); GalerkinSingleLayer / HelmholtzKernel instead drop
            % it, so the two differ by the known Delta term in the header.
            value = -alpha;
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
