function delta = bdfDelta(zeta, method)
%BDFDELTA BDF generating function delta(zeta) for Lubich convolution quadrature.
%
%   delta = bdfDelta(zeta, "BDF1" | "BDF2").  The CQ Laplace nodes are
%   s = bdfDelta(zeta)/dt with zeta on the rho-circle.  BDF1 = 1 - zeta;
%   BDF2 = 3/2 - 2 zeta + 1/2 zeta^2.  Both are A-stable, so the image of the
%   |zeta| < 1 circle lies in the right half-plane Re(s) > 0 (see the CQ solvers
%   volTdBem / volFemBemCoupled / volFemBemElastic and
%   visualizeConvolutionQuadrature).
switch upper(method)
    case "BDF1"
        delta = 1 - zeta;
    case "BDF2"
        delta = 1.5 - 2*zeta + 0.5*zeta.^2;
end
end
