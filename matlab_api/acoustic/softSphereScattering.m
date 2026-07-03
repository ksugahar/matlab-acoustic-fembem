function ref = softSphereScattering(waveNumber, radius, points, options)
%SOFTSPHERESCATTERING Partial-wave series for plane-wave sound-soft scattering.
%
%   ref = softSphereScattering(k, R, points);
%   ref.scattered      % p_scat at the points (exterior, r >= R)
%   ref.incident       % p_inc = exp(1i*k*z)
%   ref.total          % p_inc + p_scat (zero on the sphere surface)
%
% Plane wave exp(1i*k*z) hitting a sound-soft (p = 0) sphere of radius R
% centered at the origin:
%
%   p_scat(r,th) = -sum_l i^l (2l+1) [j_l(kR)/h_l(kR)] h_l(kr) P_l(cos th)
%
% with spherical Bessel j_l / Hankel h_l = h_l^(1) functions (e^{+ikr}
% radiation convention) and Legendre polynomials P_l by upward recurrence.
% The series is truncated after Terms modes (default ceil(kR) + 12, ample
% below kR ~ 10); the magnitude of the last added mode is reported as
% truncationTail so gates can verify the tail is negligible.

arguments
    waveNumber (1,1) double {mustBePositive}
    radius (1,1) double {mustBePositive}
    points (:,3) double
    options.Terms (1,1) double {mustBeInteger, mustBePositive} = ceil(waveNumber * radius) + 12
end

r = sqrt(sum(points.^2, 2));
if any(r < radius * (1 - 1e-9))
    error("softSphereScattering:interior", ...
        "Evaluation points must lie on or outside the sphere r >= R.");
end
costh = points(:, 3) ./ r;
k = waveNumber;

scattered = zeros(size(r));
Pprev = ones(size(costh));       % P_0
Pcurr = costh;                   % P_1
lastMode = zeros(size(r));
for l = 0:options.Terms
    if l == 0
        Pl = Pprev;
    elseif l == 1
        Pl = Pcurr;
    else
        Pl = ((2 * l - 1) .* costh .* Pcurr - (l - 1) .* Pprev) ./ l;
        Pprev = Pcurr;
        Pcurr = Pl;
    end
    a = -(1i^l) * (2 * l + 1) ...
        * sphericalBesselJ(l, k * radius) ./ sphericalHankel1(l, k * radius);
    lastMode = a .* sphericalHankel1(l, k * r) .* Pl;
    scattered = scattered + lastMode;
end

ref = struct();
ref.kind = "soft_sphere_plane_wave_scattering_series";
ref.policy = "analytic_partial_wave_reference_e_plus_ikr_convention";
ref.wavenumber = k;
ref.radius = radius;
ref.terms = options.Terms;
ref.truncationTail = max(abs(lastMode));
ref.scattered = scattered;
ref.incident = exp(1i * k .* points(:, 3));
ref.total = ref.incident + ref.scattered;
end


function j = sphericalBesselJ(l, x)
%SPHERICALBESSELJ j_l(x) through the half-integer cylindrical Bessel J.
j = sqrt(pi ./ (2 * x)) .* besselj(l + 0.5, x);
end


function h = sphericalHankel1(l, x)
%SPHERICALHANKEL1 h_l^(1)(x) = j_l(x) + 1i*y_l(x), outgoing for e^{+ikr}.
h = sqrt(pi ./ (2 * x)) .* (besselj(l + 0.5, x) + 1i * bessely(l + 0.5, x));
end
