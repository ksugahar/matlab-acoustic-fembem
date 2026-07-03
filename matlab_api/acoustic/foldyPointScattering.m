function res = foldyPointScattering(waveNumber, radius, centers, points)
%FOLDYPOINTSCATTERING Sound-soft spheres as coupled monopole point scatterers.
%
%   res = foldyPointScattering(k, R, centers, points);
%   res.scattered / res.incident / res.total    % at the points
%   res.excitingField                           % E_i at each center
%
% Foldy (1945) self-consistent multiple scattering with the EXACT
% single-sphere monopole (l = 0) coefficient of a sound-soft sphere:
%
%   f   = -j_0(kR) / h_0(kR)
%   E_i = p_inc(x_i) + f * sum_{j~=i} h_0(k |x_i - x_j|) E_j
%   p_scat(x) = f * sum_i E_i h_0(k |x - x_i|)
%
% Incident field: plane wave exp(1i*k*z), e^{+ikr} convention throughout.
% This is the analytic-CLASS reference for the sonic-crystal rung: closed
% form except one small N x N solve, literature-standard for sonic
% crystals. Point approximation caveats (measured in the tests, not
% assumed): only the monopole is kept, so the l >= 1 single-sphere terms
% are the truncation error, and center spacing must stay well above R.

arguments
    waveNumber (1,1) double {mustBePositive}
    radius (1,1) double {mustBePositive}
    centers (:,3) double
    points (:,3) double
end

k = waveNumber;
n = size(centers, 1);
f = -sphBesselJ0(k * radius) / sphHankel0(k * radius);

H = zeros(n, n);
for i = 1:n
    for j = 1:n
        if i ~= j
            H(i, j) = sphHankel0(k * norm(centers(i, :) - centers(j, :)));
        end
    end
end
pincAtCenters = exp(1i * k * centers(:, 3));
E = (eye(n) - f * H) \ pincAtCenters;

scattered = zeros(size(points, 1), 1);
for i = 1:n
    r = sqrt(sum((points - centers(i, :)).^2, 2));
    if any(r <= radius)
        error("foldyPointScattering:interior", ...
            "Evaluation point inside scatterer %d.", i);
    end
    scattered = scattered + f * E(i) * sphHankel0(k * r);
end

res = struct();
res.kind = "foldy_point_scatterer_soft_spheres";
res.policy = "analytic_class_monopole_multiple_scattering_reference";
res.wavenumber = k;
res.radius = radius;
res.monopoleCoefficient = f;
res.excitingField = E;
res.incident = exp(1i * k * points(:, 3));
res.scattered = scattered;
res.total = res.incident + scattered;
end


function j = sphBesselJ0(x)
%SPHBESSELJ0 j_0(x) = sin(x)/x.
j = sin(x) ./ x;
end


function h = sphHankel0(x)
%SPHHANKEL0 h_0^(1)(x) = -1i*exp(1i*x)/x, outgoing for e^{+ikr}.
h = -1i * exp(1i * x) ./ x;
end
