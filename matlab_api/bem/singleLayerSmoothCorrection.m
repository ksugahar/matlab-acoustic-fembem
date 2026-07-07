function C = singleLayerSmoothCorrection(targetPoints, sourcePoints, s, soundSpeed, sourceWeights)
%SINGLELAYERSMOOTHCORRECTION Smooth part (exp(-s r/c) - 1)/(4 pi r), weighted.
%
%   C = singleLayerSmoothCorrection(targetPoints, sourcePoints, s, c, weights)
%   C   % (nTarget x nSource): weights_j * (exp(-s r_ij/c) - 1)/(4 pi r_ij)
%
% The s-DEPENDENT half of the Laplace-domain single-layer operators
% (laplaceSingleLayerGalerkin, laplaceSingleLayerPotential).  The s-INDEPENDENT
% analytic Laplace panel matrix is assembled separately (and can be cached
% across CQ nodes, since it does not depend on s -- that caching is the CQ
% grid-radiation speedup).  Vectorized over the full target-by-source distance
% matrix; (exp(z)-1)/r uses a Taylor branch for small |z| (stable) and the finite
% limit -alpha at coincident points r = 0.

alpha = s / soundSpeed;
R = sqrt( (targetPoints(:,1) - sourcePoints(:,1).').^2 ...
        + (targetPoints(:,2) - sourcePoints(:,2).').^2 ...
        + (targetPoints(:,3) - sourcePoints(:,3).').^2 );
Z = -alpha * R;
V = (exp(Z) - 1) ./ R;                          % general (exp(z)-1)/r
small = abs(Z) < 1e-5;                          % Taylor-stable branch for small |z|
if any(small(:))
    Zs = Z(small); Rs = R(small);
    val = zeros(size(Zs)); term = ones(size(Zs));
    for k = 1:10
        term = term .* Zs / k;                  % Zs^k / k!
        val = val + term ./ Rs;                 % += (Zs^k/k!)/Rs
    end
    V(small) = val;
end
V(R == 0) = -alpha;                             % finite limit as r -> 0
C = (sourceWeights(:).' .* V) / (4 * pi);
end
