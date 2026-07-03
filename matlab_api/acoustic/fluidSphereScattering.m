function ref = fluidSphereScattering(waveNumber, radius, points, options)
%FLUIDSPHERESCATTERING Partial-wave series for a penetrable fluid sphere.
%
%   ref = fluidSphereScattering(k0, R, points, ...
%       "InteriorWavenumber", k1, "DensityRatio", rho1/rho0);
%   ref.total       % pressure at the points (interior series for r <= R,
%                   % incident + scattered outside)
%   ref.incident    % exp(1i*k0*z)
%
% Anderson (1950) transmission problem: plane wave exp(1i*k0*z) hitting a
% fluid sphere (sound speed / density contrast) centered at the origin.
% Per mode l with incident coefficient a_l = i^l (2l+1), the scattered
% A_l h_l(k0 r) and interior B_l j_l(k1 r) amplitudes solve the 2x2
% transmission match at r = R:
%
%   pressure : a_l j_l(k0 R) + A_l h_l(k0 R)              = B_l j_l(k1 R)
%   velocity : (k0/rho0) [a_l j_l'(k0 R) + A_l h_l'(k0 R)] = (k1/rho1) B_l j_l'(k1 R)
%
% (rho0 = 1; e^{+ikr} convention). For k1 = k0 and rho1 = rho0 the sphere
% is acoustically INVISIBLE: A_l = 0, B_l = a_l, total == incident - the
% exact null gate used by the coupled FEM/BEM tests. Derivatives use
% f_l'(x) = f_{l-1}(x) - (l+1)/x f_l(x), valid down to l = 0 through the
% half-integer Bessel identities (j_{-1} = cos(x)/x, h_{-1} = exp(1i*x)/x).

arguments
    waveNumber (1,1) double {mustBePositive}
    radius (1,1) double {mustBePositive}
    points (:,3) double
    options.InteriorWavenumber (1,1) double {mustBePositive} = waveNumber
    options.DensityRatio (1,1) double {mustBePositive} = 1.0
    options.Terms (1,1) double {mustBeInteger, mustBePositive} = ...
        ceil(max(waveNumber, waveNumber) * radius) + 12
end

k0 = waveNumber;
k1 = options.InteriorWavenumber;
rhor = options.DensityRatio;

r = sqrt(sum(points.^2, 2));
% the exterior partial-wave sum must converge at the FARTHEST requested
% point (k0*r terms), not just on the sphere - key for far probes
rMax = max([radius; r]);
L = max(options.Terms, ceil(max(k0 * rMax, k1 * radius)) + 12);
rSafe = max(r, 1e-30);
costh = points(:, 3) ./ rSafe;
inside = r <= radius * (1 + 1e-12);

x0 = k0 * radius;
x1 = k1 * radius;

total = zeros(size(r));
Pprev = ones(size(costh));
Pcurr = costh;
lastMode = zeros(size(r));
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

    aInc = (1i^l) * (2 * l + 1);
    % transmission match, solved by ANALYTIC elimination through the
    % interior log-derivative beta = (k1/rho1) j_l'(x1)/j_l(x1) - the naive
    % 2x2 solve is catastrophically ill-conditioned at high l (h_l huge,
    % j_l tiny) and polluted the invisible case; this form gives A = 0
    % EXACTLY when k1 = k0 and rho1 = rho0.
    j0R = sphJ(l, x0);   h0R = sphH(l, x0);   j1R = sphJ(l, x1);
    dj0R = sphJ(l - 1, x0) - (l + 1) / x0 * j0R;
    dh0R = sphH(l - 1, x0) - (l + 1) / x0 * h0R;
    dj1R = sphJ(l - 1, x1) - (l + 1) / x1 * j1R;
    beta = (k1 / rhor) * dj1R / j1R;
    A = -aInc * (k0 * dj0R - beta * j0R) / (k0 * dh0R - beta * h0R);
    B = (aInc * j0R + A * h0R) / j1R;

    mode = zeros(size(r));
    mode(inside) = B * sphJ(l, k1 * rSafe(inside)) .* Pl(inside);
    mode(~inside) = (aInc * sphJ(l, k0 * rSafe(~inside)) ...
        + A * sphH(l, k0 * rSafe(~inside))) .* Pl(~inside);
    lastMode = mode;
    total = total + mode;
end

ref = struct();
ref.kind = "fluid_sphere_transmission_scattering_series";
ref.policy = "analytic_anderson_transmission_reference_e_plus_ikr_convention";
ref.wavenumber = k0;
ref.interiorWavenumber = k1;
ref.densityRatio = rhor;
ref.radius = radius;
ref.terms = L;
ref.truncationTail = max(abs(lastMode));
ref.incident = exp(1i * k0 .* points(:, 3));
ref.total = total;
ref.insideMask = inside;
end


function j = sphJ(l, x)
%SPHJ Spherical Bessel j_l via half-integer J (j_{-1} = cos(x)/x included).
j = sqrt(pi ./ (2 * x)) .* besselj(l + 0.5, x);
end


function h = sphH(l, x)
%SPHH Spherical Hankel h_l^(1) via half-integer J/Y (h_{-1} = exp(1i*x)/x).
h = sqrt(pi ./ (2 * x)) .* (besselj(l + 0.5, x) + 1i * bessely(l + 0.5, x));
end
