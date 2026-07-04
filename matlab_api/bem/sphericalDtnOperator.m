function op = sphericalDtnOperator(surface, options)
%SPHERICALDTNOPERATOR Exact spherical Helmholtz Dirichlet-to-Neumann operator.
%
%   op = sphericalDtnOperator(surface, "Wavenumber", k);
%   op = sphericalDtnOperator(surface, "Wavenumber", k, "Degree", 8);
%   q = op.apply(p);      % nodal outward normal flux dp/dn from nodal pressure p
%   D = op.matrix;        % weak DtN <T p, mu_i> (nB x nB, rank (N+1)^2)
%
% The exact exterior Dirichlet-to-Neumann map on a SPHERE truncation.  In the
% acoustic lane this is a spherical radiating-impedance/DtN operator, not a
% Kelvin boundary.  For an outgoing solution
% p = sum a_nm h_n(kr)/h_n(kR) Y_n^m, the
% normal derivative on r = R is diagonal in spherical harmonics:
%
%   dp/dn|_Gamma = sum a_nm Lambda_n Y_n^m ,   Lambda_n = k h_n'(kR)/h_n(kR)
%
% (h_n = spherical Hankel of the first kind, e^{+ikr} convention; h_n' via the
% recurrence h_n'(x) = h_{n-1}(x) - (n+1)/x h_n(x), the SAME recurrence the
% rigid-sphere series uses). The weak DtN on the P1 boundary space is
%
%   D = Mb Phi diag(Lambda) (Phi' Mb Phi)^{-1} Phi' Mb ,   Phi_i,(n,m) = Y_n^m(dir_i)
%
% a rank-(N+1)^2 operator (N = Degree). It REPLACES the dense N^2 Galerkin
% single/double-layer assembly (no singular integration, no dense exterior
% matrix) - the fast exterior path when the truncation is a sphere. The
% general (non-spherical) radiator stays on the dense BEM: this operator is
% FAIL-LOUD and raises if the surface is not a sphere (no silent fallback).

arguments
    surface (1,1) SurfaceMesh
    options.Wavenumber (1,1) double {mustBePositive}
    options.Degree (1,1) double {mustBeInteger} = -1            % -1 = auto from kR
    options.SphericityTolerance (1,1) double {mustBePositive} = 3e-2
end

k = options.Wavenumber;
vtx = surface.vtx;
nB = size(vtx, 1);

% ---- sphericity gate via algebraic sphere fit (probe, don't guess) ----
% fit center c and radius R in one least-squares pass so the deviation
% measures true faceting, not a node-centroid offset: for each node,
%   2 x . c + (R^2 - |c|^2) = |x|^2  ->  [2 x, 1] [c; s] = |x|^2, s = R^2-|c|^2.
Afit = [2 * vtx, ones(nB, 1)];
pfit = Afit \ sum(vtx.^2, 2);
center = pfit(1:3).';
R = sqrt(pfit(4) + center * center.');
rel = vtx - center;
radii = sqrt(sum(rel.^2, 2));
deviation = max(abs(radii - R)) / R;
if deviation > options.SphericityTolerance
    error("sphericalDtnOperator:notSpherical", ...
        "DtN exterior requires a spherical truncation surface; max radius " + ...
        "deviation %.3g exceeds tol %.3g (R = %.4g, center = [%.3g %.3g %.3g]). " + ...
        "Use the dense BEM exterior for a general radiator.", ...
        deviation, options.SphericityTolerance, R, center(1), center(2), center(3));
end

% ---- degree: auto from kR, capped by the mesh resolution ((N+1)^2 <= nB) ----
N = options.Degree;
if N < 0
    N = ceil(k * R) + 8;                       % propagating n ~ kR + evanescent tail
end
maxByMesh = max(1, floor(sqrt(nB)) - 1);
if (N + 1)^2 > nB
    N = min(N, maxByMesh);
end

% ---- real spherical harmonics at the node directions ----
dirs = rel ./ radii;                           % unit directions on the sphere
[Phi, degreeOf] = realSphericalHarmonics(dirs, N);   % nB x (N+1)^2, degree tag

% ---- DtN eigenvalues Lambda_n = k h_n'(kR)/h_n(kR) ----
x0 = k * R;
Lambda = zeros(N + 1, 1);
for n = 0:N
    hn = sphH(n, x0);
    dhn = sphH(n - 1, x0) - (n + 1) / x0 * sphH(n, x0);
    Lambda(n + 1) = k * dhn / hn;
end
lamCol = Lambda(degreeOf + 1);                 % per-column eigenvalue (nModes x 1)

% ---- weak DtN  D = Mb Phi diag(Lambda) Gram^{-1} Phi' Mb ----
[Mb, ~] = SurfaceP1Space(surface).mass();
MbPhi = Mb * Phi;                              % nB x nModes
Gram = Phi' * MbPhi;                           % nModes x nModes (real SPD)
coreT = lamCol .* (Gram \ MbPhi.');            % diag(Lambda) Gram^{-1} Phi' Mb, nModes x nB
D = MbPhi * coreT;                             % nB x nB, complex-symmetric, low rank

op = struct();
op.kind = "spherical_helmholtz_dirichlet_to_neumann_operator";
op.policy = "exact_sphere_dtn_kelvin_operator_low_rank_exterior";
op.wavenumber = k;
op.radius = R;
op.center = center;
op.degree = N;
op.numModes = (N + 1)^2;
op.sphericityDeviation = deviation;
op.eigenvalues = Lambda;
op.gramCondition = cond(Gram);
op.matrix = D;
op.surfaceMass = Mb;
op.harmonics = Phi;                            % nB x nModes  Y_n^m at the nodes
op.gram = Gram;                                % nModes x nModes  Phi' Mb Phi
op.modeEigenvalues = lamCol;                   % nModes x 1  Lambda_n per column
op.apply = @(p) Phi * (lamCol .* (Gram \ (MbPhi.' * p)));   % nodal T[p] = dp/dn
op.checks = struct( ...
    "sphericalSurface", deviation <= options.SphericityTolerance, ...
    "gramWellPosed", isfinite(op.gramCondition) && op.gramCondition < 1e12, ...
    "modesResolved", (N + 1)^2 <= nB);
if all(structfun(@(v) logical(v), op.checks))
    op.status = "ok";
else
    op.status = "needs_attention";
end
end


function [Phi, degreeOf] = realSphericalHarmonics(dirs, N)
%REALSPHERICALHARMONICS Orthonormal real spherical harmonics at unit directions.
% Columns ordered n = 0..N, within n the order m = 0, then (+m cos, -m sin)
% for m = 1..n. degreeOf(col) = n so each column can be tagged with Lambda_n.

np = size(dirs, 1);
ct = dirs(:, 3);                               % cos(theta) = z on the unit sphere
ph = atan2(dirs(:, 2), dirs(:, 1));
nModes = (N + 1)^2;
Phi = zeros(np, nModes);
degreeOf = zeros(nModes, 1);
col = 1;
for n = 0:N
    P = legendre(n, ct);                       % (n+1) x np, rows m = 0..n (CS phase)
    if n == 0
        P = reshape(P, 1, np);                 % legendre(0,.) collapses to a vector
    end
    for m = 0:n
        % c_nm = sqrt((2n+1)/(4pi) * (n-m)!/(n+m)!), factorial ratio via gammaln
        cnm = sqrt((2 * n + 1) / (4 * pi)) * ...
            exp(0.5 * (gammaln(n - m + 1) - gammaln(n + m + 1)));
        Pm = P(m + 1, :).';                    % np x 1
        if m == 0
            Phi(:, col) = cnm * Pm;
            degreeOf(col) = n;
            col = col + 1;
        else
            s2 = sqrt(2) * cnm;
            Phi(:, col)     = s2 * Pm .* cos(m * ph);
            degreeOf(col)   = n;
            Phi(:, col + 1) = s2 * Pm .* sin(m * ph);
            degreeOf(col + 1) = n;
            col = col + 2;
        end
    end
end
end


function j = sphJ(l, x) %#ok<DEFNU>
%SPHJ Spherical Bessel j_l via half-integer J.
j = sqrt(pi ./ (2 * x)) .* besselj(l + 0.5, x);
end


function h = sphH(l, x)
%SPHH Spherical Hankel h_l^(1) via half-integer J/Y (h_{-1} = e^{ix}/x).
h = sqrt(pi ./ (2 * x)) .* (besselj(l + 0.5, x) + 1i * bessely(l + 0.5, x));
end
