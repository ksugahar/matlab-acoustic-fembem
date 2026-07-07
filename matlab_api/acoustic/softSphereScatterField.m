function field = softSphereScatterField(volFile, options)
%SOFTSPHERESCATTERFIELD Time-domain plane-wave pulse scattering off a soft sphere.
%
%   field = softSphereScatterField() sweeps an incident plane-wave Ricker pulse
%   in +z onto a SOFT sphere (pressure p = 0 on the surface), radiates the
%   scattered wave with the Lubich CQ time-domain BEM, and samples the total
%   field p_inc + p_scattered on an x-z plane grid at every CQ time step.  The
%   result feeds writeSoftSphereScatterGif for an animation.
%
%   The soft boundary condition sets the scattered Dirichlet trace to -incident
%   on the surface (BoundaryTimeData = -p_inc(boundary, t)).  A back-scatter
%   listener records the retroreflected arrival so the causal onset can be
%   checked against the geometric (ray) travel time.
%
%   field.pressure is [nx, nz, nt] with NaN inside the sphere; field.checks are
%   the causal/real/bounded/arrival gates; field.summary holds the scalars.

arguments
    volFile (1,1) string = "S:/MATLAB/Gypsilab/fixtures/mesh_topology/unit_sphere_coarse.vol"
    options.Radius (1,1) double {mustBePositive} = 1.0
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 3)} = 22
    options.TimeStep (1,1) double {mustBePositive} = 0.32
    options.PulseStart (1,1) double = -3.0
    options.PulseWidth (1,1) double {mustBePositive} = 0.7
    options.GridExtent (1,1) double {mustBePositive} = 4.0
    options.NumGrid (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumGrid, 3)} = 70
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder,[1 3 7])} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
end

c0 = options.SoundSpeed;
R = options.Radius;
N = options.NumTime;
dt = options.TimeStep;
t = (0:N-1).' * dt;

% --- incident plane-wave Ricker pulse travelling in +z ------------------------ %
z0 = options.PulseStart;
pw = options.PulseWidth;
ricker = @(tau) (1 - 2*(tau/pw).^2) .* exp(-(tau/pw).^2);
pInc = @(P, tt) ricker(tt - (P(:,3) - z0)/c0);            % p_inc(P, t)

% --- SOFT sphere: scattered Dirichlet trace = -incident on the boundary ------- %
mesh = VolMesh(volFile);
surface = mesh.boundary();
Xb = surface.vtx;
nB = size(Xb, 1);
boundaryData = zeros(N, nB);
for k = 1:N
    boundaryData(k, :) = -pInc(Xb, t(k)).';
end

% --- x-z plane grid of listeners (y = 0), skip points inside the sphere ------- %
ext = options.GridExtent;
ax = linspace(-ext, ext, options.NumGrid);
az = linspace(-ext, ext, options.NumGrid);
[GX, GZ] = meshgrid(ax, az);
inside = (GX(:).^2 + GZ(:).^2) < (R*1.05)^2;
gridOut = [GX(~inside), zeros(nnz(~inside),1), GZ(~inside)];

% --- back-scatter (ray) listener: +z pulse hits -z pole, retroreflects down --- %
zPole = -R;
listener = [0 0 z0];                                     % below the sphere
obsAll = [gridOut; listener];

% --- radiate the scattered field with the CQ solver --------------------------- %
cq = volTdBemConvolutionQuadrature(volFile, ...
    NumTime=N, TimeStep=dt, SoundSpeed=c0, Method=options.Method, ...
    QuadratureOrder=options.QuadratureOrder, ...
    BoundaryTimeData=boundaryData, ObservationPoints=obsAll);
scatOut = real(cq.pressure(:, 1:end-1));                 % (N x nGridOut) scattered
scatL = real(cq.pressure(:, end));                       % back-scatter listener

% --- compose total = incident + scattered on the full grid (NaN inside) ------- %
pressure = nan(options.NumGrid, options.NumGrid, N);
for k = 1:N
    total = nan(numel(GX), 1);
    total(~inside) = pInc(gridOut, t(k)) + scatOut(k, :).';
    pressure(:, :, k) = reshape(total, size(GX));
end

% --- geometric (ray) arrival + measured onset/peak of the scattered listener -- %
tHit = (zPole - z0) / c0;                                % pulse peak reaches -z pole
geoArrival = tHit + norm(listener - [0 0 zPole]) / c0;   % + retroreflect travel
peakLevel = max(abs(scatL));
onsetIdx = find(abs(scatL) > 0.1*peakLevel, 1, "first");
if isempty(onsetIdx), onsetIdx = N; end
measuredOnset = t(onsetIdx);
[~, iPk] = max(abs(scatL));
measuredPeak = t(iPk);

scale = max(1, max(abs(pressure(~isnan(pressure))), [], "all"));
field = struct();
field.kind = "soft_sphere_scatter_field";
field.volFile = string(volFile);
field.radius = R;
field.sound_speed = c0;
field.time = t;
field.time_step = dt;
field.x = ax;
field.z = az;
field.mask_inside = reshape(inside, size(GX));
field.pressure = pressure;
field.scattered_listener = scatL;
field.listener = listener;
% The scattered back-reflection PEAKS within a few steps of the geometric (ray)
% arrival: a soft sphere is not a clean point retroreflector -- the whole
% illuminated cap contributes, so the peak lags the specular first-arrival by
% ~1 step -- but it must not precede the ray (causality) nor lag it absurdly.
field.checks = struct( ...
    "finite_pressure", all(isfinite(pressure(~isnan(pressure)))), ...
    "real_time_response", cq.summary.max_imag_pressure_before_real < 1e-8*scale, ...
    "arrival_peak_causal", measuredPeak >= geoArrival - dt, ...
    "arrival_peak_near_geometry", measuredPeak <= geoArrival + 3*dt, ...
    "cq_residuals_small", cq.summary.max_relative_residual < 1e-6);
field.summary = struct( ...
    "num_time", N, ...
    "num_grid_points", size(gridOut,1), ...
    "max_abs_pressure", max(abs(pressure(~isnan(pressure))), [], "all"), ...
    "max_condition_number", cq.summary.max_condition_number, ...
    "max_relative_residual", cq.summary.max_relative_residual, ...
    "geometric_arrival", geoArrival, ...
    "measured_onset", measuredOnset, ...
    "measured_peak", measuredPeak);
if all(structfun(@(v) logical(v), field.checks))
    field.status = "ok";
else
    field.status = "needs_attention";
end
end
