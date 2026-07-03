function res = acousticFocusAdjoint(surface, sources, target, wavenumber, amplitudes, options)
%ACOUSTICFOCUSADJOINT Adjoint sensitivity of focused intensity to array phases.
%
%   res = acousticFocusAdjoint(surface, sources, target, k, amplitudes);
%   res.objective       % J = |u(target)|^2, the focused intensity
%   res.gradientConj    % dJ/d conj(p): the steepest-ascent direction
%   res.gradientReal    % dJ/d Re(p_j),  res.gradientImag = dJ/d Im(p_j)
%   res.adjointSolves   % == 1, independent of the number of sources
%
% Wavefront synthesis by REVERSE-MODE (adjoint) automatic differentiation
% through the rigid-scattering BEM solve. A phased array of point sources
% (complex amplitudes p, one per source) insonifies a sound-hard scatterer;
% the objective is the focused intensity J = |u(target)|^2 at a target point
% (e.g. focusing acoustic energy behind the scatterer - the seed of acoustic
% radiation-force / thrust design).
%
% The total field is affine in the amplitudes:
%
%   A t = M (S p),   A = 1/2 M - K_k    (rigid total-field BEM)
%   u(target) = S0 p + D_k[t](target) = (S0 + d0 A^{-1} M S) p = w p
%
% with S(n,j) = G_k(vtx_n, src_j), S0(j) = G_k(target, src_j), and
% d0 = doubleLayerPotentialMatrix(target). The row w is assembled with ONE
% adjoint (transpose) solve that is INDEPENDENT of the number of sources:
%
%   A.' lambda = d0.'  ->  w = S0 + lambda.' M S.
%
% That is the whole point of the adjoint: N design variables, one extra
% solve, versus N forward re-solves for finite differences. J is
% non-holomorphic (|.|^2), so the real gradient uses Wirtinger calculus:
%   dJ/dRe(p_j) =  2 Re(conj(u) w_j),  dJ/dIm(p_j) = -2 Im(conj(u) w_j).
% The complex steepest-ASCENT direction of a real objective is the p-bar
% Wirtinger derivative dJ/dconj(p) = 2 u conj(w) (equal, elementwise, to
% gradientReal + 1i*gradientImag) - NOT 2 conj(u) w (that is dJ/dp, the
% conjugate / wrong-way direction; near a zero field almost any step raises
% |u|^2 so the sign only bites away from the minimum).
%
% With GradientCheck=true the same gradient is verified against central
% finite differences (the repo's standard gradient-check discipline; the
% forward map is exact-linear so the check lands at ~1e-10).

arguments
    surface (1,1) SurfaceMesh
    sources (:,3) double
    target (1,3) double
    wavenumber (1,1) double {mustBePositive}
    amplitudes (:,1) double
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 7
    options.GradientCheck (1,1) logical = false
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
end

nSrc = size(sources, 1);
if numel(amplitudes) ~= nSrc
    error("acousticFocusAdjoint:amplitudes", ...
        "amplitudes must have one entry per source (%d).", nSrc);
end
p = amplitudes(:);
k = wavenumber;
vtx = surface.vtx;

% incident-field maps: surface trace and target value as linear in p
Ssurf = pointSourceMatrix(vtx, sources, k);      % nNodes x nSrc
S0 = pointSourceMatrix(target, sources, k);      % 1 x nSrc

[M, ~] = SurfaceP1Space(surface).mass();
K = GalerkinDoubleLayer(surface, "Wavenumber", k, ...
    "QuadratureOrder", options.QuadratureOrder).matrix;
A = 0.5 * M - K;
d0 = doubleLayerPotentialMatrix(surface, target, k, options.QuadratureOrder);

% forward solve (rigid total-field trace), then the affine field value
[Lf, Uf, Pf] = lu(A);
t = Uf \ (Lf \ (Pf * (M * (Ssurf * p))));
u = S0 * p + d0 * t;

% adjoint solve: ONE transpose solve, independent of nSrc
lambda = A.' \ d0.';
w = S0 + (lambda.' * M) * Ssurf;                 % 1 x nSrc, u == w p exactly

objective = abs(u)^2;
gradientReal = 2 * real(conj(u) * w).';          % dJ/dRe(p_j), N x 1
gradientImag = -2 * imag(conj(u) * w).';         % dJ/dIm(p_j), N x 1
ascentDirection = gradientReal + 1i * gradientImag;   % = 2 u conj(w), N x 1

res = struct();
res.kind = "acoustic_phased_array_focus_adjoint_sensitivity";
res.policy = "reverse_mode_adjoint_through_rigid_bem_solve";
res.wavenumber = k;
res.target = target;
res.numSources = nSrc;
res.field = u;
res.objective = objective;
res.sensitivityRow = w;
res.gradientReal = gradientReal;
res.gradientImag = gradientImag;
res.ascentDirection = ascentDirection;
res.forwardLinearityResidual = abs(u - w * p);
res.adjointSolves = 1;
res.quadratureOrder = options.QuadratureOrder;

checks = struct("forwardAffineExact", res.forwardLinearityResidual < 1e-10);

if options.GradientCheck
    h = options.FiniteDifferenceStep;
    Jof = @(pp) abs(S0 * pp + d0 * (Uf \ (Lf \ (Pf * (M * (Ssurf * pp))))))^2;
    fdRe = zeros(nSrc, 1);
    fdIm = zeros(nSrc, 1);
    for j = 1:nSrc
        e = zeros(nSrc, 1);
        e(j) = 1;
        fdRe(j) = (Jof(p + h * e) - Jof(p - h * e)) / (2 * h);
        fdIm(j) = (Jof(p + 1i * h * e) - Jof(p - 1i * h * e)) / (2 * h);
    end
    res.finiteDifferenceReal = fdRe;
    res.finiteDifferenceImag = fdIm;
    res.gradientCheckRelError = max( ...
        norm(gradientReal - fdRe) / max(norm(fdRe), eps), ...
        norm(gradientImag - fdIm) / max(norm(fdIm), eps));
    checks.gradientMatchesFiniteDifference = res.gradientCheckRelError < 1e-6;
end

res.checks = checks;
if all(structfun(@(v) logical(v), checks))
    res.status = "ok";
else
    res.status = "needs_attention";
end
end


function G = pointSourceMatrix(points, sources, k)
%POINTSOURCEMATRIX G(i,j) = exp(1i k |x_i - s_j|)/(4 pi |x_i - s_j|).
nP = size(points, 1);
nS = size(sources, 1);
G = zeros(nP, nS);
for j = 1:nS
    r = sqrt(sum((points - sources(j, :)).^2, 2));
    if any(r <= eps)
        error("acousticFocusAdjoint:coincident", ...
            "A point coincides with source %d.", j);
    end
    G(:, j) = exp(1i * k .* r) ./ (4 * pi * r);
end
end
