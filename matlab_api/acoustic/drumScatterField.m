function field = drumScatterField(volFile, options)
%DRUMSCATTERFIELD Time-domain x-z field of a two-spot drum roll scattered by a sphere.
%
%   field = drumScatterField() strikes the DRUMHEAD (top face) of a CYLINDRICAL
%   drum ALTERNATELY at two spots (+x on odd beats, -x on even), radiates the
%   drum roll with the Lubich CQ time-domain BEM, and lets the radiated sound be
%   SCATTERED by a SPHERE floating directly above the drum RIM.  The drum shell,
%   drumhead, and scatterer are carried by ONE combined .vol and ONE CQ single-
%   layer solve, so the sampled field contains the drum radiation AND the
%   multiple scattering off the sphere (both bodies enforce the sound-soft
%   total-pressure boundary).  The field is sampled on an x-z plane grid at every
%   CQ time step.
%
%   This is the MATLAB-ONLY (no Gmsh) movie companion of the drum + scatterer
%   scene: it returns a field struct shaped for writeSoftSphereScatterGif, so the
%   whole animation is produced inside MATLAB (indexed-image GIF, headless).
%
%   The two bodies are auto-detected from the combined surface (connected
%   components of the boundary triangulation); the DRUM is the body reaching
%   lowest in z, the SCATTERER is the body floating above.  No geometry is hard-
%   coded -- the drum radius, drumhead height, and scatterer center/radius are all
%   measured from the mesh.  field.pressure is [nz, nx, nt] with NaN inside either
%   body's cross-section.
%
%   FRAME COUNT / CQ RADIUS: identical to drumRollField -- the renderer uses a
%   movie-stable radius rho = eps^(0.3/2N) (rho^-N ~ 2e2, N-independent) so the
%   tail stays physical at large frame counts.  Pass CqRadius to override.

arguments
    volFile (1,1) string = "S:/MATLAB/Gypsilab/fixtures/mesh_topology/drum_scatter.vol"
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.NumBeats (1,1) double {mustBeInteger, mustBePositive} = 3
    options.BeatInterval (1,1) double {mustBePositive} = 1.0
    options.TimeStep (1,1) double {mustBePositive} = 0.06
    options.TailTime (1,1) double {mustBePositive} = 4.0
    options.StrikeWidthTime (1,1) double {mustBePositive} = 0.35
    options.StrikeWidthSpace (1,1) double {mustBePositive} = 0.33
    options.StrikeOffset (1,1) double {mustBePositive} = 0.6
    options.GridExtent (1,1) double {mustBePositive} = 1.75
    options.NumGrid (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumGrid, 3)} = 90
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder,[1 3 7])} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
    options.CqRadius double = []
end

c0 = options.SoundSpeed;

% --- combined surface, split into the two bodies (drum + scatterer) ------------ %
surface = VolMesh(volFile).boundary();
X = surface.vtx; nB = size(X, 1);
comp = surfaceConnectedComponents(nB, surface.tri);   % node -> body label
numBodies = max(comp);
loZ = accumarray(comp, X(:, 3), [numBodies 1], @min);
[~, drumBody] = min(loZ);                              % drum reaches lowest
drumMask = comp == drumBody;
scatMask = ~drumMask;

zTop = max(X(drumMask, 3)); zBot = min(X(drumMask, 3));
headTol = 0.05 * (zTop - zBot + eps);
topMask = drumMask & (abs(X(:, 3) - zTop) < max(headTol, 1e-6));
Rd = max(vecnorm(X(topMask, 1:2), 2, 2));             % drum radius (measured)
scCenter = mean(X(scatMask, :), 1);                   % scatterer center (measured)
scR = max(vecnorm(X(scatMask, :) - scCenter, 2, 2));  % scatterer radius (measured)

% --- two strike spots on the drumhead (top face) at +x and -x ----------------- %
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

% --- alternating localized Ricker taps: A on odd beats, B on even beats -------- %
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

% --- x-z plane grid (y=0) covering both bodies; skip both cross-sections ------- %
ext = options.GridExtent;
xlo = -Rd - ext; xhi = max(Rd, scCenter(1) + scR) + ext;
zlo = zBot - 0.75*ext; zhi = max(zTop, scCenter(3) + scR) + ext;
ax = linspace(xlo, xhi, options.NumGrid);
az = linspace(zlo, zhi, options.NumGrid);
[GX, GZ] = meshgrid(ax, az);
inDrum = (abs(GX(:)) <= Rd*1.06) & (GZ(:) >= zBot - 0.05) & (GZ(:) <= zTop + 0.05);
inScat = (GX(:) - scCenter(1)).^2 + (GZ(:) - scCenter(3)).^2 <= (scR*1.06)^2;
masked = inDrum | inScat;
gridOut = [GX(~masked), zeros(nnz(~masked),1), GZ(~masked)];

% --- ONE CQ solve on the combined surface (strong coupling; multiple scatter) -- %
cq = volTdBemConvolutionQuadrature(volFile, ...
    NumTime=N, TimeStep=dt, SoundSpeed=c0, Method=options.Method, ...
    QuadratureOrder=options.QuadratureOrder, CqRadius=cqRadius, ...
    BoundaryTimeData=boundaryData, ObservationPoints=gridOut);
scatOut = real(cq.pressure);

pressure = nan(options.NumGrid, options.NumGrid, N);
for k = 1:N
    frame = nan(numel(GX), 1);
    frame(~masked) = scatOut(k, :).';
    pressure(:, :, k) = reshape(frame, size(GX));
end

scale = max(1, max(abs(pressure(~isnan(pressure))), [], "all"));
field = struct();
field.kind = "drum_scatter_two_body_field";
field.volFile = string(volFile);
field.drum_radius = Rd;
field.scatterer_center = scCenter;
field.scatterer_radius = scR;
field.sound_speed = c0;
field.time = t;
field.time_step = dt;
field.x = ax;
field.z = az;
field.mask_inside = reshape(masked, size(GX));
field.pressure = pressure;
field.strikeSpotA = spotA;
field.strikeSpotB = spotB;
field.beatSpot = beatSpot;
field.checks = struct( ...
    "finite_pressure", all(isfinite(pressure(~isnan(pressure)))), ...
    "real_time_response", cq.summary.max_imag_pressure_before_real < 1e-8*scale, ...
    "cq_residuals_small", cq.summary.max_relative_residual < 1e-6, ...
    "alternating_two_spot", numel(unique(beatSpot)) == 2 && beatSpot(1) == "A", ...
    "two_bodies_detected", numBodies == 2, ...
    "scatterer_above_drum", min(X(scatMask, 3)) > zTop);
field.summary = struct( ...
    "num_time", N, ...
    "num_frames", N, ...
    "num_grid_points", size(gridOut, 1), ...
    "num_beats", options.NumBeats, ...
    "cq_radius", cqRadius, ...
    "max_abs_pressure", max(abs(pressure(~isnan(pressure))), [], "all"), ...
    "max_condition_number", cq.summary.max_condition_number, ...
    "max_relative_residual", cq.summary.max_relative_residual, ...
    "drum_radius", Rd, ...
    "scatterer_center", scCenter, ...
    "scatterer_radius", scR);
if all(structfun(@(v) logical(v), field.checks))
    field.status = "ok";
else
    field.status = "needs_attention";
end
end


function comp = surfaceConnectedComponents(nNodes, tri)
% Label each node 1..K by connected component of the boundary triangulation.
% Two disjoint closed surfaces (drum, scatterer) -> two components.
E = sort([tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])], 2);
E = unique(E, "rows");
G = graph(E(:, 1), E(:, 2), [], nNodes);
comp = conncomp(G).';
end


function tap = rickerTap(t, centerTime, width)
x = (t - centerTime) / width;
tap = (1 - 2*x.^2) .* exp(-x.^2);
tap(abs(tap) < 1e-14) = 0;
end
