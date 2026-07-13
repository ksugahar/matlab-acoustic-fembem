function field = drumRollField(volFile, options)
%DRUMROLLFIELD Time-domain x-z field of a two-spot alternating cylinder drum roll.
%
%   field = drumRollField() strikes the DRUMHEAD (the top face) of a CYLINDRICAL
%   drum ALTERNATELY at two spots (+x on odd beats, -x on even), radiates the
%   drum roll with the Lubich CQ time-domain BEM, and samples the radiated
%   pressure on an x-z plane grid at every CQ time step.  It is the SPATIAL-FIELD
%   (movie) companion of the listener-time-series drumRollConvolutionQuadrature,
%   shaped for writeSoftSphereScatterGif.
%
%   The drum sits with its axis along z (bottom z=0, drumhead z=H); the wavefronts
%   radiate up and around from the alternating +x/-x drumhead strikes.
%   field.pressure is [nz, nx, nt] with NaN inside the drum cross-section.
%
%   FRAME COUNT / CQ RADIUS: the default TimeStep = 0.03 gives ~200 frames; the
%   renderer uses a movie-stable radius rho = eps^(0.3/2N) (rho^-N ~ 2e2,
%   N-independent) so the tail stays physical at large frame counts (~0.5%
%   aliasing hidden by the auto-scaled colormap).  Pass CqRadius to override.

arguments
    volFile (1,1) string = ""
    options.Radius (1,1) double {mustBePositive} = 1.0
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.NumBeats (1,1) double {mustBeInteger, mustBePositive} = 3
    options.BeatInterval (1,1) double {mustBePositive} = 1.0
    options.TimeStep (1,1) double {mustBePositive} = 0.03
    options.TailTime (1,1) double {mustBePositive} = 3.0
    options.StrikeWidthTime (1,1) double {mustBePositive} = 0.35
    options.StrikeWidthSpace (1,1) double {mustBePositive} = 0.33
    options.StrikeOffset (1,1) double {mustBePositive} = 0.6
    options.GridExtent (1,1) double {mustBePositive} = 3.5
    options.NumGrid (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumGrid, 3)} = 90
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder,[1 3 7])} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
    options.CqRadius double = []
end

if strlength(volFile) == 0
    volFile = defaultFixture("drum_cylinder.vol");
end

c0 = options.SoundSpeed;
R = options.Radius;

mesh = VolMesh(volFile);
surface = mesh.boundary();
X = surface.vtx; nB = size(X, 1);
zTop = max(X(:, 3));
zBot = min(X(:, 3));

% --- two strike spots on the drumhead (top face) at +x and -x ----------------- %
headTol = 0.05 * (zTop - zBot + eps);
topMask = abs(X(:, 3) - zTop) < max(headTol, 1e-6);
spotA = [ options.StrikeOffset 0 zTop];
spotB = [-options.StrikeOffset 0 zTop];
wA = exp(-sum((X - spotA).^2, 2) / (2*options.StrikeWidthSpace^2)) .* topMask;
wB = exp(-sum((X - spotB).^2, 2) / (2*options.StrikeWidthSpace^2)) .* topMask;

dt = options.TimeStep;
finalTime = options.NumBeats*options.BeatInterval + options.TailTime;
N = max(8, ceil(finalTime / dt));
t = (0:N-1).' * dt;

if isempty(options.CqRadius)
    cqRadius = eps^(0.3 / (2 * N));
else
    cqRadius = options.CqRadius;
end

% --- alternating localized Ricker taps: A on odd beats, B on even beats ------- %
boundaryData = zeros(N, nB);
beatSpot = strings(options.NumBeats, 1);
for b = 1:options.NumBeats
    tap = rickerTap(t, b * options.BeatInterval, options.StrikeWidthTime);
    if mod(b, 2) == 1
        boundaryData = boundaryData + tap .* wA.'; beatSpot(b) = "A";
    else
        boundaryData = boundaryData + tap .* wB.'; beatSpot(b) = "B";
    end
end

% --- x-z plane grid (y=0); skip nodes inside the drum cross-section ----------- %
ext = options.GridExtent;
ax = linspace(-ext, ext, options.NumGrid);
az = linspace(zBot - 1.2, zTop + ext - 0.4, options.NumGrid);
[GX, GZ] = meshgrid(ax, az);
inDrum = (abs(GX(:)) <= R*1.06) & (GZ(:) >= zBot - 0.05) & (GZ(:) <= zTop + 0.05);
gridOut = [GX(~inDrum), zeros(nnz(~inDrum),1), GZ(~inDrum)];

% --- radiate the drum roll with the (mass-consistent) CQ solver --------------- %
cq = volTdBemConvolutionQuadrature(volFile, ...
    NumTime=N, TimeStep=dt, SoundSpeed=c0, Method=options.Method, ...
    QuadratureOrder=options.QuadratureOrder, CqRadius=cqRadius, ...
    BoundaryTimeData=boundaryData, ObservationPoints=gridOut);
scatOut = real(cq.pressure);

pressure = nan(options.NumGrid, options.NumGrid, N);
for k = 1:N
    frame = nan(numel(GX), 1);
    frame(~inDrum) = scatOut(k, :).';
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
field.mask_inside = reshape(inDrum, size(GX));
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


function volFile = defaultFixture(name)
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end
