function field = drumRollField(volFile, options)
%DRUMROLLFIELD Time-domain x-z field of a two-spot alternating drum roll (CQ).
%
%   field = drumRollField() strikes the sphere ALTERNATELY at two spots (the +x
%   and -x poles), radiates the drum roll with the Lubich CQ time-domain BEM, and
%   samples the radiated pressure on an x-z plane grid at every CQ time step.  It
%   is the SPATIAL-FIELD (movie) companion of the listener-time-series
%   drumRollConvolutionQuadrature, shaped for writeSoftSphereScatterGif.
%
%   Spot A fires on beats 1,3,5,..., spot B on 2,4,6,..., so the wavefronts
%   radiate alternately from the +x then the -x pole.  field.pressure is
%   [nz, nx, nt] with NaN inside the sphere; field.checks are the finite/real/
%   residual/alternating gates.
%
%   FRAME COUNT / CQ RADIUS.  The number of movie frames IS the number of CQ
%   time steps N = ceil(finalTime / TimeStep); the default TimeStep = 0.03 gives
%   ~200 frames (smooth playback AND enough time resolution to resolve the
%   StrikeWidthTime taps and the moving wavefronts).  At the CQ solver's
%   theoretically-optimal radius rho = eps^(1/(2N)) the rho^-n unscaling
%   amplifies round-off by 1/sqrt(eps) ~ 7e7 at the last step, which blows the
%   late-time tail up for large N.  A movie only needs ~visual accuracy, so this
%   renderer uses a slightly LARGER, N-INDEPENDENT radius rho = eps^(0.3/(2N))
%   by default: rho^-N = eps^-0.15 ~ 2e2 (stable to arbitrarily many frames) at
%   the cost of ~rho^N = eps^0.15 ~ 5e-3 aliasing -- invisible under the
%   auto-scaled colormap.  Pass an explicit CqRadius in (0,1) to override.

arguments
    volFile (1,1) string = "S:/MATLAB/Gypsilab/fixtures/mesh_topology/unit_sphere_coarse.vol"
    options.Radius (1,1) double {mustBePositive} = 1.0
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.NumBeats (1,1) double {mustBeInteger, mustBePositive} = 3
    options.BeatInterval (1,1) double {mustBePositive} = 1.0
    options.TimeStep (1,1) double {mustBePositive} = 0.03
    options.TailTime (1,1) double {mustBePositive} = 3.0
    options.StrikeWidthTime (1,1) double {mustBePositive} = 0.35
    options.StrikeWidthSpace (1,1) double {mustBePositive} = 0.6
    options.GridExtent (1,1) double {mustBePositive} = 3.5
    options.NumGrid (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumGrid, 3)} = 90
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder,[1 3 7])} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
    options.CqRadius double = []
end

c0 = options.SoundSpeed;
R = options.Radius;

mesh = VolMesh(volFile);
surface = mesh.boundary();
X = surface.vtx; nB = size(X, 1);

% --- two strike spots: surface nodes furthest along +x (A) and -x (B) --------- %
[~, iA] = max(X(:,1)); spotA = X(iA, :);
[~, iB] = min(X(:,1)); spotB = X(iB, :);
wA = exp(-sum((X - spotA).^2, 2) / (2*options.StrikeWidthSpace^2));
wB = exp(-sum((X - spotB).^2, 2) / (2*options.StrikeWidthSpace^2));

dt = options.TimeStep;
finalTime = options.NumBeats*options.BeatInterval + options.TailTime;
N = max(8, ceil(finalTime / dt));
t = (0:N-1).' * dt;

% movie-stable CQ radius (N-independent round-off amplification); see header
if isempty(options.CqRadius)
    cqRadius = eps^(0.3 / (2 * N));
else
    cqRadius = options.CqRadius;
end

% --- alternating localized Ricker taps: A on odd beats, B on even beats ------- %
boundaryData = zeros(N, nB);
beatSpot = strings(options.NumBeats, 1);
for b = 1:options.NumBeats
    tb = b * options.BeatInterval;
    tap = rickerTap(t, tb, options.StrikeWidthTime);
    if mod(b, 2) == 1
        boundaryData = boundaryData + tap .* wA.'; beatSpot(b) = "A";
    else
        boundaryData = boundaryData + tap .* wB.'; beatSpot(b) = "B";
    end
end

% --- x-z plane grid (poles on the horizontal axis); skip inside the sphere ---- %
ext = options.GridExtent;
ax = linspace(-ext, ext, options.NumGrid);
az = linspace(-ext, ext, options.NumGrid);
[GX, GZ] = meshgrid(ax, az);
inside = (GX(:).^2 + GZ(:).^2) < (R*1.05)^2;
gridOut = [GX(~inside), zeros(nnz(~inside),1), GZ(~inside)];

% --- radiate the drum roll with the (mass-consistent) CQ solver --------------- %
cq = volTdBemConvolutionQuadrature(volFile, ...
    NumTime=N, TimeStep=dt, SoundSpeed=c0, Method=options.Method, ...
    QuadratureOrder=options.QuadratureOrder, CqRadius=cqRadius, ...
    BoundaryTimeData=boundaryData, ObservationPoints=gridOut);
scatOut = real(cq.pressure);

pressure = nan(options.NumGrid, options.NumGrid, N);
for k = 1:N
    frame = nan(numel(GX), 1);
    frame(~inside) = scatOut(k, :).';
    pressure(:, :, k) = reshape(frame, size(GX));
end

scale = max(1, max(abs(pressure(~isnan(pressure))), [], "all"));
field = struct();
field.kind = "drum_roll_two_spot_field";
field.volFile = string(volFile);
field.radius = R;
field.sound_speed = c0;
field.time = t;
field.time_step = dt;
field.x = ax;
field.z = az;
field.mask_inside = reshape(inside, size(GX));
field.pressure = pressure;
field.strikeSpotA = spotA;
field.strikeSpotB = spotB;
field.beatSpot = beatSpot;
field.checks = struct( ...
    "finite_pressure", all(isfinite(pressure(~isnan(pressure)))), ...
    "real_time_response", cq.summary.max_imag_pressure_before_real < 1e-8*scale, ...
    "cq_residuals_small", cq.summary.max_relative_residual < 1e-6, ...
    "alternating_two_spot", numel(unique(beatSpot)) == 2 && beatSpot(1) == "A");
field.summary = struct( ...
    "num_time", N, ...
    "num_frames", N, ...
    "num_grid_points", size(gridOut, 1), ...
    "num_beats", options.NumBeats, ...
    "cq_radius", cqRadius, ...
    "max_abs_pressure", max(abs(pressure(~isnan(pressure))), [], "all"), ...
    "max_condition_number", cq.summary.max_condition_number, ...
    "max_relative_residual", cq.summary.max_relative_residual);
if all(structfun(@(v) logical(v), field.checks))
    field.status = "ok";
else
    field.status = "needs_attention";
end
end


function tap = rickerTap(t, centerTime, width)
x = (t - centerTime) / width;
tap = (1 - 2*x.^2) .* exp(-x.^2);
tap(abs(tap) < 1e-14) = 0;
end
