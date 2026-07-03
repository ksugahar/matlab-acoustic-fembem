function ref = rigidSphereScattering(waveNumber, radius, points, options)
%RIGIDSPHERESCATTERING Partial-wave series for plane-wave sound-hard scattering.
%
%   ref = rigidSphereScattering(k, R, points);
%   ref.scattered      % A_l h_l(kr) P_l sum at the (exterior) points
%   ref.total          % incident + scattered (dp/dn = 0 on r = R)
%
% Plane wave exp(1i*k*z) hitting a sound-hard (rigid, dp/dn = 0) sphere:
%
%   A_l = -i^l (2l+1) j_l'(kR) / h_l'(kR)
%
% with derivatives through f_l'(x) = f_{l-1}(x) - (l+1)/x f_l(x) (valid at
% l = 0 via the half-integer identities). e^{+ikr} convention. The series
% term count keys on the FARTHEST evaluation point (k * r_max), and the
% last-mode magnitude is reported as truncationTail.

arguments
    waveNumber (1,1) double {mustBePositive}
    radius (1,1) double {mustBePositive}
    points (:,3) double
    options.Terms (1,1) double {mustBeInteger, mustBePositive} = 1
end

k = waveNumber;
r = sqrt(sum(points.^2, 2));
if any(r < radius * (1 - 1e-9))
    error("rigidSphereScattering:interior", ...
        "Evaluation points must lie on or outside the sphere r >= R.");
end
costh = points(:, 3) ./ r;
L = max(options.Terms, ceil(k * max([radius; r])) + 12);

scattered = zeros(size(r));
Pprev = ones(size(costh));
Pcurr = costh;
lastMode = zeros(size(r));
x0 = k * radius;
for l = 0:L
    if l == 0
        Pl = Pprev;
    elseif l == 1
        Pl = Pcurr;
    else
        Pl = ((2 * l - 1) .* costh .* Pcurr - (l - 1) .* Pprev) ./ l;
        Pprev = Pcurr;
        Pcurr = Pl;
    end
    dj = sphJ(l - 1, x0) - (l + 1) / x0 * sphJ(l, x0);
    dh = sphH(l - 1, x0) - (l + 1) / x0 * sphH(l, x0);
    A = -(1i^l) * (2 * l + 1) * dj / dh;
    lastMode = A * sphH(l, k * r) .* Pl;
    scattered = scattered + lastMode;
end

ref = struct();
ref.kind = "rigid_sphere_plane_wave_scattering_series";
ref.policy = "analytic_partial_wave_reference_e_plus_ikr_convention";
ref.wavenumber = k;
ref.radius = radius;
ref.terms = L;
ref.truncationTail = max(abs(lastMode));
ref.incident = exp(1i * k .* points(:, 3));
ref.scattered = scattered;
ref.total = ref.incident + scattered;
end


function j = sphJ(l, x)
%SPHJ Spherical Bessel j_l via half-integer J (j_{-1} = cos(x)/x included).
j = sqrt(pi ./ (2 * x)) .* besselj(l + 0.5, x);
end


function h = sphH(l, x)
%SPHH Spherical Hankel h_l^(1) via half-integer J/Y (h_{-1} = exp(1i*x)/x).
h = sqrt(pi ./ (2 * x)) .* (besselj(l + 0.5, x) + 1i * bessely(l + 0.5, x));
end
