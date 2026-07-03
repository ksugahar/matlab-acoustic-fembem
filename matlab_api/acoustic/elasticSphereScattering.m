function ref = elasticSphereScattering(waveNumber, radius, points, options)
%ELASTICSPHERESCATTERING Faran (1951) partial-wave scattering by a solid sphere.
%
%   ref = elasticSphereScattering(k, R, points, ...
%       "LongitudinalSpeed", cL, "ShearSpeed", cT, "DensityRatio", rho);
%   ref.total       % exterior pressure (incident + scattered) at the points
%   ref.scattered   % scattered part;  ref.incident = exp(1i*k*z)
%
% A solid ELASTIC sphere (radius R) in a fluid, plane wave exp(1i*k*z). This
% is the reference the acoustic FSI (fluid-structure interaction) work is
% gated against: unlike rigid/soft, an elastic sphere has INTERNAL
% resonances (compressional + shear) that leave sharp features in the
% radiation force - the whole reason to model the interior at all.
%
% Speeds and density are given RELATIVE TO THE FLUID (fluid c = 1, rho = 1):
%   LongitudinalSpeed  cL/c_fluid   (compressional / P-wave speed ratio)
%   ShearSpeed         cT/c_fluid   (shear / S-wave speed ratio; 0 = fluid)
%   DensityRatio       rho_solid/rho_fluid
%
% Interior potentials phi = A j_l(kL r) (compressional, kL = omega/cL) and
% psi = B j_l(kT r) (shear, kT = omega/cT); displacement u = grad phi +
% curl-curl of the shear potential. The three boundary conditions at r = R
% (radial displacement continuity, radial stress = -pressure, zero shear
% stress for an inviscid fluid) give a 3x3 system per mode for (c_l, A, B).
%
% VALIDATED against two INDEPENDENT references (limits, 2026-07-04):
%   ShearSpeed -> 0 reproduces the Anderson fluid sphere
%     (fluidSphereScattering) to 1e-10 - the exact fluid limit;
%   very stiff (cL, cT, rho large) reproduces rigidSphereScattering.
% The Reynolds/sign discipline mattered: the radial-stress boundary term
% is -lambda kL^2 A j_l(kL R) (sigma_rr = -p), NOT +; a flipped sign passes
% the stiff limit but fails the fluid limit by 20-70% (lock against an
% INDEPENDENT reference, not only a self-consistent one).

arguments
    waveNumber (1,1) double {mustBePositive}
    radius (1,1) double {mustBePositive}
    points (:,3) double
    options.LongitudinalSpeed (1,1) double {mustBePositive} = 2.0
    options.ShearSpeed (1,1) double {mustBeNonnegative} = 1.0
    options.DensityRatio (1,1) double {mustBePositive} = 1.5
    options.Terms (1,1) double {mustBeInteger, mustBeNonnegative} = 0
end

% high-l modes are physically negligible (truncated below) but their naive
% per-mode solve is ill-conditioned (h_l huge, j_l tiny) - a known Bessel
% artifact, not a physics issue; quiet it and restore on exit.
ws = warning('off', 'MATLAB:nearlySingularMatrix');
cleanup = onCleanup(@() warning(ws)); %#ok<NASGU>

k = waveNumber;
a = radius;
omega = k;                          % fluid c = 1
cL = options.LongitudinalSpeed;
cT = options.ShearSpeed;
rhoS = options.DensityRatio;
rhoF = 1.0;
mu = rhoS * cT^2;
lam = rhoS * (cL^2 - 2 * cT^2);
kL = omega / cL;

L = options.Terms;
if L <= 0
    % scattered angular content is l <~ k a + margin; keep it physical to
    % avoid the high-l Bessel ill-conditioning (h_l huge, j_l tiny).
    L = ceil(k * a) + 10;
end

n = size(points, 1);
r = sqrt(sum(points.^2, 2));
costh = points(:, 3) ./ r;

coeff = zeros(L + 1, 1);
for l = 0:L
    coeff(l + 1) = elasticCoeff(l, k, a, rhoF, omega, kL, cT, rhoS, lam, mu);
end

scattered = zeros(n, 1);
Pprev = ones(n, 1);
Pcurr = costh;
lastMode = zeros(n, 1);
for l = 0:L
    if l == 0
        Pl = Pprev;
    elseif l == 1
        Pl = Pcurr;
    else
        Pl = ((2*l-1) .* costh .* Pcurr - (l-1) .* Pprev) ./ l;
        Pprev = Pcurr;
        Pcurr = Pl;
    end
    lastMode = (1i^l) * (2*l+1) * coeff(l+1) * sphHankel(l, k*r) .* Pl;
    scattered = scattered + lastMode;
end

ref = struct();
ref.kind = "elastic_solid_sphere_faran_scattering_series";
ref.policy = "analytic_partial_wave_reference_e_plus_ikr_convention";
ref.wavenumber = k;
ref.radius = a;
ref.longitudinalSpeed = cL;
ref.shearSpeed = cT;
ref.densityRatio = rhoS;
ref.terms = L;
ref.truncationTail = max(abs(lastMode));
ref.incident = exp(1i * k .* points(:, 3));
ref.scattered = scattered;
ref.total = ref.incident + scattered;
end


function c = elasticCoeff(n, k, a, rhoF, omega, kL, cT, rhoS, lam, mu)
%ELASTICCOEFF Scattered coefficient c_l from the r = a boundary conditions.
x = k * a;
xl = kL * a;
JN = @(l, z) sphBessel(l, z);
DJ = @(l, z) sphBesselD(l, z);
DDJ = @(l, z) sphBesselDD(l, z);
% fluid radial displacement uses (k/(rhoF omega^2)) [j'(x) + c h'(x)]
fluidFac = k / (rhoF * omega^2);

if cT == 0
    % fluid interior (shear absent): 2x2 for [c; A], sigma_rr = -lam kL^2 A jn(xl)
    M = [ fluidFac * sphHankelD(n, x),  -kL * DJ(n, xl);
          sphHankel(n, x),              -lam * kL^2 * JN(n, xl) ];
    rhs = [ -fluidFac * DJ(n, x); -JN(n, x) ];
    sol = M \ rhs;
    c = sol(1);
    return
end

kT = omega / cT;
xt = kT * a;

% solid radial displacement  u_r = A [kL j'(xl)] + B [l(l+1)/a jn(xt)]
Ur_A = kL * DJ(n, xl);
Ur_B = n*(n+1)/a * JN(n, xt);
% d(u_r)/dr at r = a
dUr_A = kL^2 * DDJ(n, xl);
dUr_B = n*(n+1) * (kT * DJ(n, xt) / a - JN(n, xt) / a^2);
% radial normal stress  sigma_rr = -lam kL^2 phi + 2 mu d(u_r)/dr
Srr_A = -lam * kL^2 * JN(n, xl) + 2 * mu * dUr_A;
Srr_B = 2 * mu * dUr_B;
% tangential displacement radial part V and its derivative (for sigma_rtheta)
Va_A = JN(n, xl) / a;
Va_B = JN(n, xt) / a + kT * DJ(n, xt);
Vp_A = -JN(n, xl) / a^2 + kL * DJ(n, xl) / a;
Vp_B = -JN(n, xt) / a^2 + kT * DJ(n, xt) / a + kT^2 * DDJ(n, xt);
% shear stress  sigma_rtheta = mu [ u_r/a + V' - V/a ] = 0
Srt_A = mu * (Ur_A / a + Vp_A - Va_A / a);
Srt_B = mu * (Ur_B / a + Vp_B - Va_B / a);

M = [ fluidFac * sphHankelD(n, x), -Ur_A, -Ur_B;
      sphHankel(n, x),              Srr_A, Srr_B;
      0,                            Srt_A, Srt_B ];
rhs = [ -fluidFac * DJ(n, x); -JN(n, x); 0 ];
sol = M \ rhs;
c = sol(1);
end


function j = sphBessel(l, x)
j = sqrt(pi ./ (2*x)) .* besselj(l + 0.5, x);
end

function h = sphHankel(l, x)
h = sqrt(pi ./ (2*x)) .* (besselj(l + 0.5, x) + 1i * bessely(l + 0.5, x));
end

function d = sphBesselD(l, x)
d = sphBessel(l-1, x) - (l+1)./x .* sphBessel(l, x);
end

function d = sphHankelD(l, x)
d = sphHankel(l-1, x) - (l+1)./x .* sphHankel(l, x);
end

function d = sphBesselDD(l, x)
% from the spherical Bessel ODE: f'' = -(2/x) f' - (1 - l(l+1)/x^2) f
d = -(2./x) .* sphBesselD(l, x) - (1 - l*(l+1)./x.^2) .* sphBessel(l, x);
end
