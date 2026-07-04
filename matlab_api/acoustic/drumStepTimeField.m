function field = drumStepTimeField(options)
%DRUMSTEPTIMEFIELD Time-domain field from a step-struck baffled drum mode.
%
%   field = drumStepTimeField() computes an axisymmetric r-z pressure map
%   from a circular membrane hit by a step force.  This is a readable first
%   rung for vibro-acoustic FEM/BEM: one structural membrane mode supplies
%   normal acceleration, and the exterior field is the causal Rayleigh
%   retarded-potential integral.
%
%   The pressure is normalized by the modal force per mass.  Use
%   plotDrumStepTimeField(field, k) to visualize snapshot k.

arguments
    options.Radius (1,1) double {mustBePositive} = 0.10
    options.SoundSpeed (1,1) double {mustBePositive} = 343.0
    options.Density (1,1) double {mustBePositive} = 1.2
    options.NaturalFrequency (1,1) double {mustBePositive} = 220.0
    options.DampingRatio (1,1) double {mustBeNonnegative} = 0.03
    options.ForcePerModalMass (1,1) double = 1.0
    options.NumRadialObservation (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumRadialObservation, 1)} = 80
    options.NumAxialObservation (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumAxialObservation, 1)} = 80
    options.NumSourceRadial (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumSourceRadial, 2)} = 48
    options.NumSourceAzimuth (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumSourceAzimuth, 3)} = 72
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 1)} = 80
    options.RMax (1,1) double {mustBePositive} = 0.25
    options.ZMax (1,1) double {mustBePositive} = 0.35
    options.TMax (1,1) double {mustBePositive} = 1.6e-3
end

a = options.Radius;
c = options.SoundSpeed;
rho0 = options.Density;
omegaN = 2 * pi * options.NaturalFrequency;
zeta = min(options.DampingRatio, 0.999);
forcePerMass = options.ForcePerModalMass;

rObs = linspace(0, options.RMax, options.NumRadialObservation);
zObs = linspace(0.005 * a, options.ZMax, options.NumAxialObservation);
t = linspace(0, options.TMax, options.NumTime);
[Robs, Zobs] = ndgrid(rObs, zObs);

[srcR, srcTheta, srcWeight] = diskPolarQuadrature(a, ...
    options.NumSourceRadial, options.NumSourceAzimuth);
mode = drumMode01(srcR / a);
sourceWeight = srcWeight .* mode;

p = zeros(numel(rObs), numel(zObs), numel(t));
for m = 1:numel(srcR)
    xs = srcR(m) * cos(srcTheta(m));
    ys = srcR(m) * sin(srcTheta(m));
    dist = sqrt((Robs - xs).^2 + ys.^2 + Zobs.^2);
    delay = dist / c;
    for k = 1:numel(t)
        tau = t(k) - delay;
        p(:, :, k) = p(:, :, k) + sourceWeight(m) * modalStepAcceleration( ...
            tau, omegaN, zeta, forcePerMass) ./ dist;
    end
end

% Rayleigh baffled-piston pressure: p = rho0/(2*pi) * int(dv/dt / R dS).
p = rho0 / (2 * pi) * p;

field = struct();
field.kind = "drum_step_time_field_rayleigh";
field.radius = a;
field.sound_speed = c;
field.density = rho0;
field.natural_frequency = options.NaturalFrequency;
field.damping_ratio = options.DampingRatio;
field.force_per_modal_mass = forcePerMass;
field.r = rObs;
field.z = zObs;
field.t = t;
field.pressure = p;
field.pressure_normalization = "Pa for SI ForcePerModalMass; default is normalized";
field.mode = "axisymmetric membrane mode J0(alpha01*r/a), alpha01=2.4048255577";
field.boundary_integral = "causal Rayleigh retarded potential";
field.checks = struct();
field.checks.causal_initial_field_zero = max(abs(p(:, :, 1)), [], "all") < 1e-12;
field.checks.finite_pressure = all(isfinite(p), "all");
field.checks.nonzero_after_wave_arrival = max(abs(p), [], "all") > 0;
field.summary = struct();
field.summary.max_abs_pressure = max(abs(p), [], "all");
field.summary.first_nonzero_time = firstNonzeroTime(t, p);
end


function [r, theta, w] = diskPolarQuadrature(radius, nr, nth)
dr = radius / nr;
rCenters = ((1:nr).' - 0.5) * dr;
dtheta = 2 * pi / nth;
thetaCenters = ((1:nth) - 0.5) * dtheta;

[R, T] = ndgrid(rCenters, thetaCenters);
r = R(:);
theta = T(:);
w = r * dr * dtheta;
end


function y = drumMode01(rho)
alpha01 = 2.4048255577;
y = besselj(0, alpha01 * rho);
end


function a = modalStepAcceleration(tau, omegaN, zeta, forcePerMass)
a = zeros(size(tau));
active = tau >= 0;
if ~any(active, "all")
    return
end
wd = omegaN * sqrt(max(1 - zeta^2, eps));
ta = tau(active);
a(active) = forcePerMass * exp(-zeta * omegaN * ta) .* ( ...
    cos(wd * ta) - (zeta * omegaN / wd) * sin(wd * ta));
end


function t0 = firstNonzeroTime(t, p)
maxByTime = squeeze(max(max(abs(p), [], 1), [], 2));
idx = find(maxByTime > 1e-12, 1, "first");
if isempty(idx)
    t0 = NaN;
else
    t0 = t(idx);
end
end
