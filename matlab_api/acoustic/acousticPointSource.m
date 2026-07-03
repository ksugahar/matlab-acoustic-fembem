function p = acousticPointSource(waveNumber, sourcePoint, points)
%ACOUSTICPOINTSOURCE Free-space Helmholtz point source exp(1i*k*r)/(4*pi*r).
%
%   p = acousticPointSource(k, x0, points);
%
% The fundamental solution of the Helmholtz equation in the e^{+ikr}
% teaching convention (k = 0 degenerates to the real Laplace point source).
%
% This is the EXACT analytic gate for the exterior Dirichlet solve: take
% the boundary data g from a source placed INSIDE the closed surface; by
% uniqueness of the exterior Helmholtz problem with the Sommerfeld
% radiation condition, the solved field must reproduce this same function
% at every exterior point. No series truncation is involved, so the
% remaining error is pure discretization.

arguments
    waveNumber (1,1) double {mustBeNonnegative}
    sourcePoint (1,3) double
    points (:,3) double
end

r = sqrt(sum((points - sourcePoint).^2, 2));
if any(r <= eps)
    error("acousticPointSource:coincident", ...
        "Evaluation point coincides with the source point.");
end
p = exp(1i * waveNumber .* r) ./ (4 * pi * r);
if waveNumber == 0
    p = real(p);
end
end
