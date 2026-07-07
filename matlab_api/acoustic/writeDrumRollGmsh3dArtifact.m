function artifact = writeDrumRollGmsh3dArtifact(volFile, options)
%WRITEDRUMROLLGMSH3DARTIFACT Export the cylinder drum roll to a 3D Gmsh scene.
%
%   artifact = writeDrumRollGmsh3dArtifact() builds a Gmsh v4.1 post-processing
%   scene for the two-spot cylinder drum roll: the DRUMHEAD (top face of the
%   cylinder) carries a time-stepped VECTOR displacement (View VectorType=5
%   deforms the mesh, so the drumhead visibly bumps in +z at the +x/-x strike
%   spots), and the x-z plane carries the CQ-radiated acoustic pressure as a
%   scalar.  Opening the .geo animates the deforming drum + the radiating field
%   in 3D; the .msh is the raw mesh/data inspection entry point.
%
%   The scene mirrors writeDrumBemPlaneGmshArtifact (x-z pressure + a deforming
%   3D drum surface) but for the CQ TD-BEM two-spot drum roll on the real drum
%   shape -- a cylinder, not a sphere.

arguments
    volFile (1,1) string = "S:/MATLAB/Gypsilab/fixtures/mesh_topology/drum_cylinder.vol"
    options.OutputBase (1,1) string = ""
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.NumBeats (1,1) double {mustBeInteger, mustBePositive} = 3
    options.BeatInterval (1,1) double {mustBePositive} = 1.0
    options.TimeStep (1,1) double {mustBePositive} = 0.075
    options.TailTime (1,1) double {mustBePositive} = 3.0
    options.StrikeWidthTime (1,1) double {mustBePositive} = 0.35
    options.StrikeWidthSpace (1,1) double {mustBePositive} = 0.33
    options.StrikeOffset (1,1) double {mustBePositive} = 0.6
    options.GridExtent (1,1) double {mustBePositive} = 3.0
    options.NumGrid (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumGrid, 8)} = 90
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 1
    options.Method (1,1) string {mustBeMember(options.Method, ["BDF1", "BDF2"])} = "BDF2"
    options.DisplacementFactor (1,1) double {mustBePositive} = 0.55
    options.PressureRangeScale (1,1) double {mustBePositive} = 0.4
end

totalTimer = tic;
if strlength(options.OutputBase) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    options.OutputBase = fullfile(tempdir, "drum_roll_gmsh3d_" + stamp);
end

% ---------- drum surface + two-spot drumhead excitation ----------
surface = VolMesh(volFile).boundary();
X = surface.vtx; cylTri = surface.tri; nB = size(X, 1);
zTop = max(X(:, 3)); zBot = min(X(:, 3));
headTol = 0.05 * (zTop - zBot + eps);
topMask = abs(X(:, 3) - zTop) < max(headTol, 1e-6);
spotA = [ options.StrikeOffset 0 zTop];
spotB = [-options.StrikeOffset 0 zTop];
wA = exp(-sum((X - spotA).^2, 2) / (2*options.StrikeWidthSpace^2)) .* topMask;
wB = exp(-sum((X - spotB).^2, 2) / (2*options.StrikeWidthSpace^2)) .* topMask;

dt = options.TimeStep;
N = max(8, ceil((options.NumBeats*options.BeatInterval + options.TailTime) / dt));
t = (0:N-1).' * dt;
E = zeros(N, nB); beatSpot = strings(options.NumBeats, 1);
for b = 1:options.NumBeats
    tap = rickerTap(t, b*options.BeatInterval, options.StrikeWidthTime);
    if mod(b, 2) == 1, E = E + tap.*wA.'; beatSpot(b) = "A";
    else, E = E + tap.*wB.'; beatSpot(b) = "B"; end
end
disp3 = zeros(nB, 3, N);
disp3(:, 3, :) = reshape(E.', nB, 1, N);          % drumhead moves in +z

% ---------- x-z plane mesh (y=0), excluding the drum cross-section ----------
ext = options.GridExtent;
ax = linspace(-ext, ext, options.NumGrid);
az = linspace(zBot - 1.2, zTop + ext, options.NumGrid);
[GX, GZ] = meshgrid(ax, az);
Rd = max(vecnorm(X(topMask, 1:2), 2, 2));         % drum radius
inDrum = (abs(GX(:)) <= Rd*1.06) & (GZ(:) >= zBot - 0.05) & (GZ(:) <= zTop + 0.05);
outside = ~inDrum;
gi = reshape(1:numel(GX), options.NumGrid, options.NumGrid);
planeTriAll = zeros(0, 3);
for iz = 1:options.NumGrid-1
    for ix = 1:options.NumGrid-1
        ids = [gi(iz,ix), gi(iz,ix+1), gi(iz+1,ix+1), gi(iz+1,ix)];
        if all(outside(ids([1 2 3]))), planeTriAll(end+1,:) = ids([1 2 3]); end %#ok<AGROW>
        if all(outside(ids([1 3 4]))), planeTriAll(end+1,:) = ids([1 3 4]); end %#ok<AGROW>
    end
end
usedP = unique(planeTriAll(:));
remapP = zeros(numel(GX), 1); remapP(usedP) = 1:numel(usedP);
planeNodes = [GX(usedP), zeros(numel(usedP),1), GZ(usedP)];
planeTri = remapP(planeTriAll);
nP = size(planeNodes, 1);

% ---------- CQ radiate the drumhead strikes to the plane ----------
solveTimer = tic;
cq = volTdBemConvolutionQuadrature(volFile, NumTime=N, TimeStep=dt, ...
    SoundSpeed=options.SoundSpeed, Method=options.Method, ...
    QuadratureOrder=options.QuadratureOrder, CqRadius=eps^(0.3/(2*N)), ...
    BoundaryTimeData=E, ObservationPoints=planeNodes);
P = real(cq.pressure);
solveSeconds = toc(solveTimer);

% ---------- write the Gmsh scene ----------
outBase = options.OutputBase;
mshPath = outBase + ".msh";
geoPath = outBase + ".geo";
geoOptPath = geoPath + ".opt";
optPath = outBase + ".opt";
jsonPath = outBase + ".result.json";

writeTimer = tic;
writeDrumRollMsh(mshPath, planeNodes, planeTri, X, cylTri, t, P, disp3, nP);
displayRange = max(eps, options.PressureRangeScale * max(abs(P), [], "all"));
writeDrumRollGeo(geoPath, mshPath, options.DisplacementFactor, displayRange);
writeDrumRollGeo(geoOptPath, "", options.DisplacementFactor, displayRange);
copyfile(geoOptPath, optPath, "f");
writeSeconds = toc(writeTimer);

artifact = struct();
artifact.schema = "matlab-acoustic-fembem.drum-roll-gmsh3d.v1";
artifact.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
artifact.matlab_version = string(version);
artifact.vol_file = string(volFile);
artifact.gmsh_msh = mshPath;
artifact.gmsh_geo = geoPath;
artifact.gmsh_geo_opt = geoOptPath;
artifact.gmsh_msh_version = "4.1";
artifact.data_kind = "x-z plane CQ acoustic pressure scalar + deforming cylinder drumhead vector displacement";
artifact.drum_shape = "cylinder";
artifact.beat_spots = beatSpot(:).';
artifact.frames = N;
artifact.time_start = t(1);
artifact.time_end = t(end);
artifact.drum_surface_nodes = nB;
artifact.drumhead_nodes = nnz(topMask);
artifact.plane_nodes = nP;
artifact.strike_spot_a = spotA;
artifact.strike_spot_b = spotB;
artifact.drum_radius = Rd;
artifact.drum_top_z = zTop;
artifact.max_abs_pressure = max(abs(P), [], "all");
artifact.max_abs_displacement = max(abs(E), [], "all");
artifact.gmsh_displacement_factor = options.DisplacementFactor;
artifact.pressure_display_range = displayRange;
artifact.max_relative_residual = cq.summary.max_relative_residual;
artifact.checks = struct( ...
    "gmsh_v41", true, ...
    "cylinder_not_sphere", zTop > zBot + 1e-6 && nnz(topMask) > 0, ...
    "drumhead_two_spots_top_face", spotA(3) == zTop && spotB(3) == zTop, ...
    "pressure_xz_plane_written", nP > 0 && all(isfinite(P), "all"), ...
    "deforming_drum_written", max(abs(E), [], "all") > 0, ...
    "alternating_two_spot", numel(unique(beatSpot)) == 2 && beatSpot(1) == "A", ...
    "cq_residual_small", cq.summary.max_relative_residual < 1e-6, ...
    "geo_opt_written", isfile(geoOptPath));
artifact.status = "needs_attention";
if all(structfun(@(v) logical(v), artifact.checks))
    artifact.status = "ok";
end
artifact.timing = struct( ...
    "cq_solve_seconds", solveSeconds, ...
    "gmsh_write_seconds", writeSeconds, ...
    "total_seconds", toc(totalTimer));
writeString(jsonPath, jsonencode(artifact));
end


function writeDrumRollMsh(path, planeNodes, planeTri, drumNodes, drumTri, t, P, disp3, nP)
nB = size(drumNodes, 1); off = nP; N = numel(t);
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0, error("writeDrumRollGmsh3dArtifact:file", "Could not open %s", path); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");
fprintf(fid, "$PhysicalNames\n2\n2 1 ""acoustic_pressure_xz_plane""\n2 2 ""drum_cylinder_surface""\n$EndPhysicalNames\n");
fprintf(fid, "$Entities\n0 0 2 0\n");
fprintf(fid, "1 %.17g %.17g %.17g %.17g %.17g %.17g 1 1 0\n", [min(planeNodes,[],1), max(planeNodes,[],1)]);
fprintf(fid, "2 %.17g %.17g %.17g %.17g %.17g %.17g 1 2 0\n", [min(drumNodes,[],1), max(drumNodes,[],1)]);
fprintf(fid, "$EndEntities\n");
fprintf(fid, "$Nodes\n2 %d 1 %d\n", nP+nB, nP+nB);
fprintf(fid, "2 1 0 %d\n", nP); fprintf(fid, "%d\n", (1:nP).'); fprintf(fid, "%.17g %.17g %.17g\n", planeNodes.');
fprintf(fid, "2 2 0 %d\n", nB); fprintf(fid, "%d\n", (off+(1:nB)).'); fprintf(fid, "%.17g %.17g %.17g\n", drumNodes.');
fprintf(fid, "$EndNodes\n");
fprintf(fid, "$Elements\n2 %d 1 %d\n", size(planeTri,1)+size(drumTri,1), size(planeTri,1)+size(drumTri,1));
fprintf(fid, "2 1 2 %d\n", size(planeTri,1)); et = 1;
for e = 1:size(planeTri,1), fprintf(fid, "%d %d %d %d\n", et, planeTri(e,1), planeTri(e,2), planeTri(e,3)); et = et+1; end
fprintf(fid, "2 2 2 %d\n", size(drumTri,1));
for e = 1:size(drumTri,1), fprintf(fid, "%d %d %d %d\n", et, off+drumTri(e,1), off+drumTri(e,2), off+drumTri(e,3)); et = et+1; end
fprintf(fid, "$EndElements\n");
for k = 1:N
    fprintf(fid, "$NodeData\n1\n""acoustic_pressure_xz_plane""\n1\n%.17g\n3\n%d\n1\n%d\n", t(k), k-1, nP);
    fprintf(fid, "%d %.17g\n", [(1:nP); P(k,:)]);
    fprintf(fid, "$EndNodeData\n");
end
for k = 1:N
    fprintf(fid, "$NodeData\n1\n""drumhead_displacement""\n1\n%.17g\n3\n%d\n3\n%d\n", t(k), k-1, nB);
    fprintf(fid, "%d %.17g %.17g %.17g\n", ...
        [(off+(1:nB)); squeeze(disp3(:,1,k)).'; squeeze(disp3(:,2,k)).'; squeeze(disp3(:,3,k)).']);
    fprintf(fid, "$EndNodeData\n");
end
clear cleanup
end


function writeDrumRollGeo(path, mshPath, dispFactor, displayRange)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0, error("writeDrumRollGmsh3dArtifact:file", "Could not open %s", path); end
cleanup = onCleanup(@() fclose(fid));
if strlength(mshPath) > 0
    [~, base, ext] = fileparts(mshPath);
    fprintf(fid, "// Cylinder drum roll: x-z acoustic field + deforming 3D drumhead.\n");
    fprintf(fid, "Merge ""%s%s"";\n", base, ext);
end
fprintf(fid, "General.InitialModule = 5;\n");
fprintf(fid, "General.Trackball = 0;\n");
fprintf(fid, "General.RotationX = -62;\n");
fprintf(fid, "General.RotationY = 0;\n");
fprintf(fid, "General.RotationZ = -28;\n");
fprintf(fid, "General.Color.Background = {255,255,255};\n");
fprintf(fid, "General.Color.Foreground = {0,0,0};\n");
fprintf(fid, "General.Color.Text = {0,0,0};\n");
fprintf(fid, "Mesh.SurfaceFaces = 0;\n");
fprintf(fid, "Mesh.SurfaceEdges = 0;\n");
fprintf(fid, "View[0].Name = ""acoustic pressure on x-z plane"";\n");
fprintf(fid, "View[0].IntervalsType = 3;\n");
fprintf(fid, "View[0].NbIso = 40;\n");
fprintf(fid, "View[0].ColormapNumber = 2;\n");
fprintf(fid, "View[0].RangeType = 2;\n");
fprintf(fid, "View[0].CustomMin = %.17g;\n", -displayRange);
fprintf(fid, "View[0].CustomMax = %.17g;\n", displayRange);
fprintf(fid, "View[0].ShowElement = 0;\n");
fprintf(fid, "View[0].Light = 0;\n");
fprintf(fid, "View[1].Name = ""deforming drumhead"";\n");
fprintf(fid, "View[1].VectorType = 5;\n");
fprintf(fid, "View[1].DisplacementFactor = %.17g;\n", dispFactor);
fprintf(fid, "View[1].ShowElement = 1;\n");
fprintf(fid, "View[1].Light = 1;\n");
fprintf(fid, "View[1].ColormapNumber = 2;\n");
fprintf(fid, "View[1].ShowScale = 0;\n");
fprintf(fid, "PostProcessing.AnimationDelay = 0.06;\n");
fprintf(fid, "PostProcessing.AnimationCycle = 0;\n");
clear cleanup
end


function tap = rickerTap(t, centerTime, width)
x = (t - centerTime) / width;
tap = (1 - 2*x.^2) .* exp(-x.^2);
tap(abs(tap) < 1e-14) = 0;
end


function writeString(path, text)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0, error("writeDrumRollGmsh3dArtifact:file", "Could not open %s", path); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", char(text));
clear cleanup
end
