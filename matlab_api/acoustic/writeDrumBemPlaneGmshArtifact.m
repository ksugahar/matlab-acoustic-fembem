function artifact = writeDrumBemPlaneGmshArtifact(source, options)
%WRITEDRUMBEMPLANEGMSHARTIFACT Export drum BEM teaching scene to Gmsh.
%
%   artifact = writeDrumBemPlaneGmshArtifact("scene.mat", ...
%       "OutputBase", "out/drum_bem") writes a post-processing artifact where
%   the acoustic field is a BEM-style pressure distribution on the x-z plane,
%   while the drum BEM/radiation surfaces are 3-D and time-deformed by a vector
%   displacement view.  The .geo is the post display recipe; the .msh is the
%   raw mesh/data inspection entry point.

arguments
    source (1,1) string = ""
    options.OutputBase (1,1) string = ""
    options.PlaneStride (1,1) double {mustBeInteger, mustBeGreaterThan(options.PlaneStride, 0)} = 2
    options.DrumRadial (1,1) double {mustBeInteger, mustBeGreaterThan(options.DrumRadial, 2)} = 18
    options.DrumAzimuth (1,1) double {mustBeInteger, mustBeGreaterThan(options.DrumAzimuth, 7)} = 72
    options.DrumAxial (1,1) double {mustBeInteger, mustBeGreaterThan(options.DrumAxial, 2)} = 12
    options.InitialTimeFraction (1,1) double {mustBeGreaterThanOrEqual(options.InitialTimeFraction, 0), mustBeLessThanOrEqual(options.InitialTimeFraction, 1)} = 0.28
    options.PressureRangeScale (1,1) double {mustBePositive} = 0.7
    options.DeformedSurfaceAmplitude (1,1) double {mustBePositive} = 0.035
end

totalTimer = tic;
loadTimer = tic;
[scene, sourcePath] = loadScene(source);
loadSeconds = toc(loadTimer);

if strlength(options.OutputBase) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    options.OutputBase = fullfile(tempdir, "drum_bem_plane_" + stamp);
end

buildTimer = tic;
[nodes, blocks, pressureNodeTags, pressureValues, drumNodeTags, displacementValues, meta] = ...
    buildDrumBemGmshData(scene, options);
buildSeconds = toc(buildTimer);

outBase = options.OutputBase;
mshPath = outBase + ".msh";
geoPath = outBase + ".geo";
geoOptPath = geoPath + ".opt";
optPath = outBase + ".opt";
mshOptPath = mshPath + ".opt";
jsonPath = outBase + ".result.json";

writeTimer = tic;
writeMsh41(mshPath, nodes, blocks, scene.t, pressureNodeTags, pressureValues, ...
    drumNodeTags, displacementValues);
mshSeconds = toc(writeTimer);

displayRange = max(eps, options.PressureRangeScale * max(abs(pressureValues), [], "all"));
maxDisp = max(sqrt(sum(displacementValues.^2, 2)), [], "all");
if ~(isfinite(maxDisp) && maxDisp > 0)
    displacementFactor = 1.0;
else
    displacementFactor = options.DeformedSurfaceAmplitude / maxDisp;
end
initialStep = max(0, min(numel(scene.t) - 1, round(options.InitialTimeFraction * (numel(scene.t) - 1))));

geoTimer = tic;
writeGeo(geoPath, mshPath, initialStep, displayRange, displacementFactor);
writeGeoOptions(geoOptPath, initialStep, displayRange, displacementFactor);
copyfile(geoOptPath, optPath, "f");
writeMeshOptions(mshOptPath);
geoSeconds = toc(geoTimer);

jsonTimer = tic;
artifact = struct();
artifact.schema = "cae-ai-lab.drum-bem-plane-gmsh.v1";
artifact.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
artifact.matlab_version = string(version);
artifact.source_mat = sourcePath;
artifact.gmsh_msh = mshPath;
artifact.gmsh_geo = geoPath;
artifact.gmsh_geo_opt = geoOptPath;
artifact.gmsh_msh_opt = mshOptPath;
artifact.gmsh_msh_version = "4.1";
artifact.data_kind = "BEM acoustic pressure on x-z plane plus 3D deforming drum BEM surfaces";
artifact.acoustic_view = "x_z_plane_scalar_pressure";
artifact.drum_view = "top_bottom_side_bem_surfaces_vector_displacement";
artifact.frames = numel(scene.t);
artifact.time_start = scene.t(1);
artifact.time_end = scene.t(end);
artifact.nodes = size(nodes, 1);
artifact.triangles = sum(arrayfun(@(b) size(b.tri, 1), blocks));
artifact.pressure_plane_nodes = numel(pressureNodeTags);
artifact.drum_surface_nodes = numel(drumNodeTags);
artifact.drum_surface_triangles = meta.drumTriangles;
artifact.plane_stride = options.PlaneStride;
artifact.max_abs_pressure = max(abs(pressureValues), [], "all");
artifact.max_abs_physical_displacement = maxDisp;
artifact.gmsh_displacement_factor = displacementFactor;
artifact.pressure_display_range = displayRange;
artifact.initial_time_step = initialStep;
artifact.initial_time = scene.t(initialStep + 1);
artifact.force_peak_time_step = NaN;
artifact.force_peak_time = NaN;
if isfield(scene, "motion") && isfield(scene.motion, "force") && ~isempty(scene.motion.force)
    [~, forcePeakIndex] = max(abs(scene.motion.force));
    artifact.force_peak_time_step = forcePeakIndex - 1;
    artifact.force_peak_time = scene.t(forcePeakIndex);
end
artifact.checks = struct( ...
    "gmsh_v41", true, ...
    "geo_opt_exact_autoload_written", isfile(geoOptPath), ...
    "msh_opt_exact_autoload_written", isfile(mshOptPath), ...
    "pressure_xz_plane_written", numel(pressureNodeTags) > 0, ...
    "drum_3d_surface_written", numel(drumNodeTags) > 0 && meta.drumTriangles > 0, ...
    "deforming_surface_written", maxDisp > 0, ...
    "finite_pressure", all(isfinite(pressureValues), "all"), ...
    "finite_displacement", all(isfinite(displacementValues), "all"), ...
    "source_scene_ok", isfield(scene, "status") && string(scene.status) == "ok");
artifact.pass = all(structfun(@(v) logical(v), artifact.checks));
artifact.status = string("needs_attention");
if artifact.pass
    artifact.status = "ok";
end
artifact.timing = struct( ...
    "load_seconds", loadSeconds, ...
    "mesh_build_seconds", buildSeconds, ...
    "msh_write_seconds", mshSeconds, ...
    "geo_opt_write_seconds", geoSeconds, ...
    "json_write_seconds", toc(jsonTimer), ...
    "total_seconds", toc(totalTimer));
writeString(jsonPath, jsonencode(artifact));
end


function [scene, sourcePath] = loadScene(source)
if strlength(source) == 0
    scene = drumFemBemCoupledDemo();
    sourcePath = "generated:drumFemBemCoupledDemo";
    return
end
S = load(source);
if ~isfield(S, "scene")
    error("writeDrumBemPlaneGmshArtifact:source", ...
        "MAT file must contain a variable named scene: %s", source);
end
scene = S.scene;
sourcePath = source;
end


function [nodes, blocks, pressureTags, pressureValues, drumTags, dispValues, meta] = ...
        buildDrumBemGmshData(scene, options)
validateScene(scene);

[planeNodes, planeTri, pressureValues] = buildPressurePlane(scene, options.PlaneStride);
[topNodes, topTri, topMeta] = buildDisk(scene.geometry.radius, 0.0, ...
    options.DrumRadial, options.DrumAzimuth, "top");
[bottomNodes, bottomTri, bottomMeta] = buildDisk(scene.geometry.radius, -scene.geometry.depth, ...
    options.DrumRadial, options.DrumAzimuth, "bottom");
[sideNodes, sideTri, sideMeta] = buildSide(scene.geometry.outer_radius, scene.geometry.depth, ...
    options.DrumAxial, options.DrumAzimuth);

blocks = struct([]);
[nodes, blocks] = appendBlock([], blocks, 1, 1, "bem_pressure_xz_plane", planeNodes, planeTri);
[nodes, blocks] = appendBlock(nodes, blocks, 2, 2, "drum_top_head_bem_surface", topNodes, topTri);
[nodes, blocks] = appendBlock(nodes, blocks, 3, 3, "drum_bottom_head_bem_surface", bottomNodes, bottomTri);
[nodes, blocks] = appendBlock(nodes, blocks, 4, 4, "drum_side_shell_bem_surface", sideNodes, sideTri);

pressureTags = blocks(1).nodeTags(:);
drumTags = [blocks(2).nodeTags(:); blocks(3).nodeTags(:); blocks(4).nodeTags(:)];
dispValues = buildDisplacement(scene, topMeta, bottomMeta, sideMeta, ...
    blocks(2).nodeTags, blocks(3).nodeTags, blocks(4).nodeTags);

meta = struct();
meta.drumTriangles = size(topTri, 1) + size(bottomTri, 1) + size(sideTri, 1);
end


function validateScene(scene)
required = ["x", "z", "t", "pressure", "masks", "geometry", "motion"];
for name = required
    if ~isfield(scene, name)
        error("writeDrumBemPlaneGmshArtifact:scene", "scene.%s is required.", name);
    end
end
end


function [nodes, tri, values] = buildPressurePlane(scene, stride)
xIdx = 1:stride:numel(scene.x);
zIdx = 1:stride:numel(scene.z);
x = scene.x(xIdx);
z = scene.z(zIdx);
[X, Z] = meshgrid(x, z);
allNodes = [X(:), zeros(numel(X), 1), Z(:)];
mask = scene.masks.boundary_domain(zIdx, xIdx) ...
    & ~scene.masks.drum_frame(zIdx, xIdx) ...
    & ~scene.masks.interior_air(zIdx, xIdx);
grid = reshape(1:numel(X), size(X));
triAll = zeros(0, 3);
for iz = 1:numel(z)-1
    for ix = 1:numel(x)-1
        ids = [grid(iz, ix), grid(iz, ix+1), grid(iz+1, ix+1), grid(iz+1, ix)];
        if all(mask(ids([1 2 3])))
            triAll(end + 1, :) = ids([1 2 3]); %#ok<AGROW>
        end
        if all(mask(ids([1 3 4])))
            triAll(end + 1, :) = ids([1 3 4]); %#ok<AGROW>
        end
    end
end
used = unique(triAll(:));
remap = zeros(size(allNodes, 1), 1);
remap(used) = 1:numel(used);
nodes = allNodes(used, :);
tri = remap(triAll);

nt = numel(scene.t);
values = zeros(numel(used), nt);
for k = 1:nt
    frame = scene.pressure(zIdx, xIdx, k);
    values(:, k) = frame(used);
end
end


function [nodes, tri, meta] = buildDisk(radius, z, nr, nth, kind)
nodes = zeros(1 + nr * nth, 3);
meta.r = zeros(size(nodes, 1), 1);
meta.theta = zeros(size(nodes, 1), 1);
nodes(1, :) = [0 0 z];
idx = @(ir, it) 1 + (ir - 1) * nth + it;
for ir = 1:nr
    r = radius * ir / nr;
    for it = 1:nth
        th = 2*pi*(it - 1) / nth;
        id = idx(ir, it);
        nodes(id, :) = [r*cos(th), r*sin(th), z];
        meta.r(id) = r;
        meta.theta(id) = th;
    end
end
tri = zeros(0, 3);
for it = 1:nth
    jt = mod(it, nth) + 1;
    if kind == "top"
        tri(end + 1, :) = [1, idx(1, it), idx(1, jt)]; %#ok<AGROW>
    else
        tri(end + 1, :) = [1, idx(1, jt), idx(1, it)]; %#ok<AGROW>
    end
end
for ir = 2:nr
    for it = 1:nth
        jt = mod(it, nth) + 1;
        a = idx(ir-1, it); b = idx(ir-1, jt);
        c = idx(ir, it); d = idx(ir, jt);
        if kind == "top"
            tri(end + 1, :) = [a, c, d]; %#ok<AGROW>
            tri(end + 1, :) = [a, d, b]; %#ok<AGROW>
        else
            tri(end + 1, :) = [a, d, c]; %#ok<AGROW>
            tri(end + 1, :) = [a, b, d]; %#ok<AGROW>
        end
    end
end
meta.kind = kind;
meta.radius = radius;
end


function [nodes, tri, meta] = buildSide(radius, depth, nz, nth)
nodes = zeros((nz + 1) * nth, 3);
meta.theta = zeros(size(nodes, 1), 1);
meta.axial = zeros(size(nodes, 1), 1);
idx = @(iz, it) iz * nth + it;
for iz = 0:nz
    s = iz / nz;
    z = -depth + s * depth;
    for it = 1:nth
        th = 2*pi*(it - 1) / nth;
        id = idx(iz, it);
        nodes(id, :) = [radius*cos(th), radius*sin(th), z];
        meta.theta(id) = th;
        meta.axial(id) = s;
    end
end
tri = zeros(0, 3);
for iz = 0:nz-1
    for it = 1:nth
        jt = mod(it, nth) + 1;
        a = idx(iz, it); b = idx(iz, jt);
        c = idx(iz + 1, it); d = idx(iz + 1, jt);
        tri(end + 1, :) = [a, c, d]; %#ok<AGROW>
        tri(end + 1, :) = [a, d, b]; %#ok<AGROW>
    end
end
meta.kind = "side";
meta.radius = radius;
end


function values = buildDisplacement(scene, topMeta, bottomMeta, sideMeta, topTags, bottomTags, sideTags)
nt = numel(scene.t);
n = numel(topTags) + numel(bottomTags) + numel(sideTags);
values = zeros(n, 3, nt);
topMode = drumMode01(topMeta.r / topMeta.radius);
bottomMode = drumMode01(bottomMeta.r / bottomMeta.radius);
sideMode = sin(pi * sideMeta.axial);
rowTop = 1:numel(topTags);
rowBottom = numel(topTags) + (1:numel(bottomTags));
rowSide = numel(topTags) + numel(bottomTags) + (1:numel(sideTags));
for k = 1:nt
    values(rowTop, 3, k) = -topMode * scene.motion.top_displacement(k);
    values(rowBottom, 3, k) = -bottomMode * scene.motion.bottom_displacement(k);
    radial = scene.motion.shell_displacement(k) * sideMode;
    values(rowSide, 1, k) = radial .* cos(sideMeta.theta);
    values(rowSide, 2, k) = radial .* sin(sideMeta.theta);
end
end


function y = drumMode01(rho)
alpha01 = 2.4048255577;
y = besselj(0, alpha01 * rho);
end


function [allNodes, blocks] = appendBlock(allNodes, blocks, entityTag, physTag, name, nodes, tri)
offset = size(allNodes, 1);
nodeTags = offset + (1:size(nodes, 1));
block = struct();
block.entityTag = entityTag;
block.physTag = physTag;
block.name = string(name);
block.nodes = nodes;
block.tri = tri + offset;
block.nodeTags = nodeTags;
if isempty(blocks)
    blocks = block;
else
    blocks(end + 1) = block; %#ok<AGROW>
end
allNodes = [allNodes; nodes]; %#ok<AGROW>
end


function writeMsh41(path, nodes, blocks, time, pressureTags, pressureValues, drumTags, dispValues)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeDrumBemPlaneGmshArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");
fprintf(fid, "$PhysicalNames\n%d\n", numel(blocks));
for b = 1:numel(blocks)
    fprintf(fid, "2 %d ""%s""\n", blocks(b).physTag, blocks(b).name);
end
fprintf(fid, "$EndPhysicalNames\n");
fprintf(fid, "$Entities\n0 0 %d 0\n", numel(blocks));
for b = 1:numel(blocks)
    bb = [min(blocks(b).nodes, [], 1), max(blocks(b).nodes, [], 1)];
    fprintf(fid, "%d %.17g %.17g %.17g %.17g %.17g %.17g 1 %d 0\n", ...
        blocks(b).entityTag, bb(1), bb(2), bb(3), bb(4), bb(5), bb(6), blocks(b).physTag);
end
fprintf(fid, "$EndEntities\n");

fprintf(fid, "$Nodes\n%d %d 1 %d\n", numel(blocks), size(nodes, 1), size(nodes, 1));
for b = 1:numel(blocks)
    fprintf(fid, "2 %d 0 %d\n", blocks(b).entityTag, numel(blocks(b).nodeTags));
    fprintf(fid, "%d\n", blocks(b).nodeTags);
    for j = 1:size(blocks(b).nodes, 1)
        fprintf(fid, "%.17g %.17g %.17g\n", blocks(b).nodes(j, 1), blocks(b).nodes(j, 2), blocks(b).nodes(j, 3));
    end
end
fprintf(fid, "$EndNodes\n");

numElements = sum(arrayfun(@(b) size(b.tri, 1), blocks));
fprintf(fid, "$Elements\n%d %d 1 %d\n", numel(blocks), numElements, numElements);
elemTag = 1;
for b = 1:numel(blocks)
    fprintf(fid, "2 %d 2 %d\n", blocks(b).entityTag, size(blocks(b).tri, 1));
    for e = 1:size(blocks(b).tri, 1)
        fprintf(fid, "%d %d %d %d\n", elemTag, blocks(b).tri(e, 1), blocks(b).tri(e, 2), blocks(b).tri(e, 3));
        elemTag = elemTag + 1;
    end
end
fprintf(fid, "$EndElements\n");

writeScalarNodeData(fid, "bem_pressure_xz_plane", time, pressureTags, pressureValues);
writeVectorNodeData(fid, "drum_surface_displacement", time, drumTags, dispValues);
clear cleanup
end


function writeScalarNodeData(fid, name, time, tags, values)
for k = 1:numel(time)
    fprintf(fid, "$NodeData\n1\n""%s""\n1\n%.17g\n3\n%d\n1\n%d\n", ...
        name, time(k), k - 1, numel(tags));
    for i = 1:numel(tags)
        fprintf(fid, "%d %.17g\n", tags(i), values(i, k));
    end
    fprintf(fid, "$EndNodeData\n");
end
end


function writeVectorNodeData(fid, name, time, tags, values)
for k = 1:numel(time)
    fprintf(fid, "$NodeData\n1\n""%s""\n1\n%.17g\n3\n%d\n3\n%d\n", ...
        name, time(k), k - 1, numel(tags));
    for i = 1:numel(tags)
        fprintf(fid, "%d %.17g %.17g %.17g\n", tags(i), values(i, 1, k), values(i, 2, k), values(i, 3, k));
    end
    fprintf(fid, "$EndNodeData\n");
end
end


function writeGeo(path, mshPath, initialStep, displayRange, displacementFactor)
[~, base, ext] = fileparts(mshPath);
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeDrumBemPlaneGmshArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Post view: BEM acoustic field on x-z plane + deforming 3D drum BEM surfaces.\n");
fprintf(fid, "Merge ""%s%s"";\n", base, ext);
writeGeoOptionsBody(fid, initialStep, displayRange, displacementFactor);
clear cleanup
end


function writeGeoOptions(path, initialStep, displayRange, displacementFactor)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeDrumBemPlaneGmshArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Exact .geo.opt post-display sidecar for drum BEM artifact.\n");
writeGeoOptionsBody(fid, initialStep, displayRange, displacementFactor);
clear cleanup
end


function writeGeoOptionsBody(fid, initialStep, displayRange, displacementFactor)
fprintf(fid, "General.InitialModule = 5;\n");
fprintf(fid, "General.Trackball = 0;\n");
fprintf(fid, "General.RotationX = -68;\n");
fprintf(fid, "General.RotationY = 0;\n");
fprintf(fid, "General.RotationZ = 0;\n");
fprintf(fid, "General.RotationCenterGravity = 1;\n");
fprintf(fid, "General.ScaleX = 1.25;\n");
fprintf(fid, "General.ScaleY = 1.25;\n");
fprintf(fid, "General.ScaleZ = 1.25;\n");
fprintf(fid, "General.Color.Background = {255,255,255};\n");
fprintf(fid, "General.Color.Foreground = {0,0,0};\n");
fprintf(fid, "General.Color.Text = {0,0,0};\n");
fprintf(fid, "Mesh.SurfaceFaces = 0;\n");
fprintf(fid, "Mesh.SurfaceEdges = 0;\n");
fprintf(fid, "Mesh.VolumeEdges = 0;\n");
fprintf(fid, "PostProcessing.Link = 1;\n");
fprintf(fid, "View[0].Visible = 1;\n");
fprintf(fid, "View[0].Name = ""BEM pressure on x-z plane"";\n");
fprintf(fid, "View[0].TimeStep = %d;\n", initialStep);
fprintf(fid, "View[0].IntervalsType = 2;\n");
fprintf(fid, "View[0].NbIso = 40;\n");
fprintf(fid, "View[0].ColormapNumber = 2;\n");
fprintf(fid, "View[0].RangeType = 2;\n");
fprintf(fid, "View[0].CustomMin = %.17g;\n", -displayRange);
fprintf(fid, "View[0].CustomMax = %.17g;\n", displayRange);
fprintf(fid, "View[0].ShowElement = 0;\n");
fprintf(fid, "View[0].ShowScale = 1;\n");
fprintf(fid, "View[0].Light = 0;\n");
fprintf(fid, "View[1].Visible = 1;\n");
fprintf(fid, "View[1].Name = ""3D deforming drum BEM surface"";\n");
fprintf(fid, "View[1].TimeStep = %d;\n", initialStep);
fprintf(fid, "View[1].VectorType = 5;\n");
fprintf(fid, "View[1].DisplacementFactor = %.17g;\n", displacementFactor);
fprintf(fid, "View[1].ShowElement = 1;\n");
fprintf(fid, "View[1].ShowScale = 0;\n");
fprintf(fid, "View[1].IntervalsType = 2;\n");
fprintf(fid, "View[1].ColormapNumber = 10;\n");
fprintf(fid, "View[1].DrawSkinOnly = 1;\n");
fprintf(fid, "View[1].Light = 1;\n");
fprintf(fid, "View[1].LightTwoSide = 1;\n");
fprintf(fid, "View[1].SmoothNormals = 1;\n");
fprintf(fid, "PostProcessing.AnimationDelay = 0.04;\n");
fprintf(fid, "PostProcessing.AnimationCycle = 0;\n");
fprintf(fid, "PostProcessing.AnimationStep = 1;\n");
end


function writeMeshOptions(path)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeDrumBemPlaneGmshArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Mesh/data inspection sidecar. Open .geo for the post display.\n");
fprintf(fid, "General.InitialModule = 1;\n");
fprintf(fid, "General.Trackball = 0;\n");
fprintf(fid, "General.RotationX = -68;\n");
fprintf(fid, "General.RotationY = 0;\n");
fprintf(fid, "General.RotationZ = 0;\n");
fprintf(fid, "Mesh.SurfaceFaces = 1;\n");
fprintf(fid, "Mesh.SurfaceEdges = 1;\n");
fprintf(fid, "Mesh.VolumeEdges = 0;\n");
fprintf(fid, "Mesh.ColorCarousel = 2;\n");
fprintf(fid, "View[0].Visible = 0;\n");
fprintf(fid, "View[1].Visible = 0;\n");
clear cleanup
end


function writeString(path, text)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeDrumBemPlaneGmshArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", char(text));
clear cleanup
end
