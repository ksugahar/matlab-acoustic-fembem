function artifact = writeGmshPostLaunchArtifact(mshPath, options)
%WRITEGMSHPOSTLAUNCHARTIFACT Write a .geo launch recipe and .opt sidecars.
%
%   artifact = writeGmshPostLaunchArtifact("case.msh", "Views", views)
%   writes:
%       case.geo        normal post-processing launch target
%       case.geo.opt    exact option sidecar auto-loaded by Gmsh
%       case.msh.opt    raw mesh/data inspection sidecar
%       case.opt        compatibility mirror only
%       case.display.json
%
%   The function is intentionally display-only: it does not generate a mesh
%   and does not edit the .msh data.  Use it when a Gypsilab-style readable
%   solver has already produced Gmsh MSH v4.1 NodeData/ElementData and the
%   result needs a durable section/camera/view recipe.

arguments
    mshPath (1,1) string
    options.OutputBase (1,1) string = ""
    options.Title (1,1) string = "Gmsh post launch artifact"
    options.CameraPreset (1,1) string {mustBeMember(options.CameraPreset, ...
        ["z_up_xz_from_positive_y", "positive_y_oblique", "front_xz", "custom"])} = "z_up_xz_from_positive_y"
    options.Rotation (1,3) double = [NaN NaN NaN]
    options.Scale (1,3) double {mustBePositive} = [1.25 1.25 1.25]
    options.EnableCutPlane (1,1) logical = false
    options.CutPlaneNormal (1,3) double = [0 -1 0]
    options.CutPlaneOffset (1,1) double = 0
    options.CutPlaneWholeElements (1,1) logical = false
    options.CutPlaneOnlyVolume (1,1) logical = false
    options.CutPlaneOnlyDrawIntersectingVolume (1,1) logical = false
    options.MeshSurfaceFaces (1,1) logical = false
    options.MeshSurfaceEdges (1,1) logical = false
    options.MeshVolumeFaces (1,1) logical = false
    options.MeshVolumeEdges (1,1) logical = false
    options.MeshNumSubEdges (1,1) double {mustBeInteger, mustBePositive} = 4
    options.LinkTimeSteps (1,1) logical = true
    options.AnimationDelay (1,1) double {mustBePositive} = 0.04
    options.AnimationCycle (1,1) logical = false
    options.AnimationStep (1,1) double {mustBeInteger, mustBePositive} = 1
    options.Views = struct([])
end

if strlength(mshPath) == 0
    error("writeGmshPostLaunchArtifact:mshPath", "mshPath is required.");
end

views = options.Views;
if ~isstruct(views)
    error("writeGmshPostLaunchArtifact:views", "Views must be a struct array.");
end
if options.EnableCutPlane && norm(options.CutPlaneNormal) == 0
    error("writeGmshPostLaunchArtifact:cutPlane", "CutPlaneNormal must be nonzero.");
end

[mshFolder, mshBase, mshExt] = fileparts(mshPath);
if strlength(mshExt) == 0
    mshExt = ".msh";
end
if strlength(options.OutputBase) == 0
    outBase = string(fullfile(mshFolder, mshBase));
else
    outBase = options.OutputBase;
end

geoPath = outBase + ".geo";
geoOptPath = geoPath + ".opt";
optPath = outBase + ".opt";
mshOptPath = string(fullfile(mshFolder, mshBase + mshExt + ".opt"));
jsonPath = outBase + ".display.json";

rotation = cameraRotation(options.CameraPreset, options.Rotation);
mergeTarget = gmshMergeTarget(mshPath, geoPath);
maxViewIndex = maxViewIndexFrom(views);

writeGeo(geoPath, mergeTarget, options, rotation, views);
writeGeoOptions(geoOptPath, options, rotation, views);
copyfile(geoOptPath, optPath, "f");
writeMeshOptions(mshOptPath, options, rotation, maxViewIndex);

artifact = struct();
artifact.schema = "cae-ai-lab.gmsh-post-launch.v1";
artifact.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
artifact.matlab_version = string(version);
artifact.gmsh_msh = string(mshPath);
artifact.gmsh_geo = geoPath;
artifact.gmsh_geo_opt = geoOptPath;
artifact.gmsh_opt = optPath;
artifact.gmsh_msh_opt = mshOptPath;
artifact.display_json = jsonPath;
artifact.launch_target = geoPath;
artifact.mesh_inspection_target = string(mshPath);
artifact.camera_preset = options.CameraPreset;
artifact.camera_rotation = rotation;
artifact.camera_axis_up = "z";
artifact.cut_plane_enabled = options.EnableCutPlane;
artifact.cut_plane_normal = options.CutPlaneNormal;
artifact.cut_plane_offset = options.CutPlaneOffset;
artifact.view_count = numel(views);
artifact.view_names = viewNames(views);
artifact.checks = struct( ...
    "geo_launch_written", isfile(geoPath), ...
    "geo_opt_exact_autoload_written", isfile(geoOptPath), ...
    "msh_opt_exact_autoload_written", isfile(mshOptPath), ...
    "plain_opt_is_compatibility_only", isfile(optPath), ...
    "launch_target_is_geo", endsWith(geoPath, ".geo"), ...
    "camera_has_z_up_contract", true, ...
    "cut_plane_metadata_recorded", ~options.EnableCutPlane || norm(options.CutPlaneNormal) > 0, ...
    "view_metadata_recorded", numel(views) > 0);
artifact.status = "needs_attention";
artifact.pass = all(structfun(@(x) logical(x), artifact.checks));
if artifact.pass
    artifact.status = "ok";
end

writeString(jsonPath, jsonencode(artifact));
end


function writeGeo(path, mergeTarget, options, rotation, views)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeGmshPostLaunchArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// %s\n", char(options.Title));
fprintf(fid, "// Open this .geo for post-processing; open .msh for raw mesh inspection.\n");
fprintf(fid, "Merge ""%s"";\n", escapeGmshString(mergeTarget));
writeOptionsBody(fid, options, rotation, views);
clear cleanup
end


function writeGeoOptions(path, options, rotation, views)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeGmshPostLaunchArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Exact .geo.opt post-display sidecar.\n");
writeOptionsBody(fid, options, rotation, views);
clear cleanup
end


function writeMeshOptions(path, options, rotation, maxViewIndex)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeGmshPostLaunchArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Mesh/data inspection sidecar. Open .geo for the post display.\n");
fprintf(fid, "General.InitialModule = 1;\n");
writeCamera(fid, rotation, options.Scale);
fprintf(fid, "General.Color.Background = {255,255,255};\n");
fprintf(fid, "General.Color.Foreground = {0,0,0};\n");
fprintf(fid, "General.Color.Text = {0,0,0};\n");
fprintf(fid, "Mesh.NumSubEdges = %d;\n", options.MeshNumSubEdges);
fprintf(fid, "Mesh.SurfaceFaces = 1;\n");
fprintf(fid, "Mesh.SurfaceEdges = 1;\n");
fprintf(fid, "Mesh.VolumeFaces = 0;\n");
fprintf(fid, "Mesh.VolumeEdges = 0;\n");
fprintf(fid, "Mesh.ColorCarousel = 2;\n");
for idx = 0:maxViewIndex
    fprintf(fid, "View[%d].Visible = 0;\n", idx);
end
clear cleanup
end


function writeOptionsBody(fid, options, rotation, views)
fprintf(fid, "General.InitialModule = 5;\n");
writeCamera(fid, rotation, options.Scale);
fprintf(fid, "General.Color.Background = {255,255,255};\n");
fprintf(fid, "General.Color.Foreground = {0,0,0};\n");
fprintf(fid, "General.Color.Text = {0,0,0};\n");
fprintf(fid, "Mesh.NumSubEdges = %d;\n", options.MeshNumSubEdges);
fprintf(fid, "Mesh.SurfaceFaces = %d;\n", double(options.MeshSurfaceFaces));
fprintf(fid, "Mesh.SurfaceEdges = %d;\n", double(options.MeshSurfaceEdges));
fprintf(fid, "Mesh.VolumeFaces = %d;\n", double(options.MeshVolumeFaces));
fprintf(fid, "Mesh.VolumeEdges = %d;\n", double(options.MeshVolumeEdges));
fprintf(fid, "Mesh.Light = 1;\n");
fprintf(fid, "Mesh.LightTwoSide = 1;\n");
writeCutPlane(fid, options);
fprintf(fid, "PostProcessing.Link = %d;\n", double(options.LinkTimeSteps));
for k = 1:numel(views)
    writeView(fid, views(k), k - 1, options.EnableCutPlane);
end
fprintf(fid, "PostProcessing.AnimationDelay = %.17g;\n", options.AnimationDelay);
fprintf(fid, "PostProcessing.AnimationCycle = %d;\n", double(options.AnimationCycle));
fprintf(fid, "PostProcessing.AnimationStep = %d;\n", options.AnimationStep);
end


function writeCamera(fid, rotation, scale)
fprintf(fid, "General.Trackball = 0;\n");
fprintf(fid, "General.RotationX = %.17g;\n", rotation(1));
fprintf(fid, "General.RotationY = %.17g;\n", rotation(2));
fprintf(fid, "General.RotationZ = %.17g;\n", rotation(3));
fprintf(fid, "General.RotationCenterGravity = 1;\n");
fprintf(fid, "General.ScaleX = %.17g;\n", scale(1));
fprintf(fid, "General.ScaleY = %.17g;\n", scale(2));
fprintf(fid, "General.ScaleZ = %.17g;\n", scale(3));
fprintf(fid, "General.TranslationX = 0;\n");
fprintf(fid, "General.TranslationY = 0;\n");
fprintf(fid, "General.TranslationZ = 0;\n");
end


function writeCutPlane(fid, options)
if ~options.EnableCutPlane
    fprintf(fid, "Mesh.Clip = 0;\n");
    return
end
n = options.CutPlaneNormal(:).' / norm(options.CutPlaneNormal);
fprintf(fid, "General.Clip0A = %.17g;\n", n(1));
fprintf(fid, "General.Clip0B = %.17g;\n", n(2));
fprintf(fid, "General.Clip0C = %.17g;\n", n(3));
fprintf(fid, "General.Clip0D = %.17g;\n", options.CutPlaneOffset);
fprintf(fid, "General.ClipFactor = 5;\n");
fprintf(fid, "General.ClipWholeElements = %d;\n", double(options.CutPlaneWholeElements));
fprintf(fid, "General.ClipOnlyVolume = %d;\n", double(options.CutPlaneOnlyVolume));
fprintf(fid, "General.ClipOnlyDrawIntersectingVolume = %d;\n", ...
    double(options.CutPlaneOnlyDrawIntersectingVolume));
fprintf(fid, "Mesh.Clip = 1;\n");
end


function writeView(fid, view, fallbackIndex, cutPlaneEnabled)
idx = double(getView(view, "Index", fallbackIndex));
name = string(getView(view, "Name", "view_" + string(idx)));
kind = lower(string(getView(view, "Kind", "scalar")));
fprintf(fid, "View[%d].Visible = %d;\n", idx, double(getView(view, "Visible", true)));
fprintf(fid, "View[%d].Name = ""%s"";\n", idx, escapeGmshString(name));
fprintf(fid, "View[%d].TimeStep = %d;\n", idx, double(getView(view, "TimeStep", 0)));
fprintf(fid, "View[%d].IntervalsType = %d;\n", idx, double(getView(view, "IntervalsType", 2)));
fprintf(fid, "View[%d].NbIso = %d;\n", idx, double(getView(view, "NbIso", 32)));
fprintf(fid, "View[%d].ColormapNumber = %d;\n", idx, double(getView(view, "ColormapNumber", 2)));
fprintf(fid, "View[%d].ShowElement = %d;\n", idx, double(getView(view, "ShowElement", false)));
fprintf(fid, "View[%d].ShowScale = %d;\n", idx, double(getView(view, "ShowScale", true)));
fprintf(fid, "View[%d].Light = %d;\n", idx, double(getView(view, "Light", false)));
fprintf(fid, "View[%d].Clip = %d;\n", idx, double(getView(view, "Clip", cutPlaneEnabled)));
range = getView(view, "Range", []);
if isnumeric(range) && numel(range) == 2
    fprintf(fid, "View[%d].RangeType = 2;\n", idx);
    fprintf(fid, "View[%d].CustomMin = %.17g;\n", idx, range(1));
    fprintf(fid, "View[%d].CustomMax = %.17g;\n", idx, range(2));
end
if kind == "vector"
    fprintf(fid, "View[%d].VectorType = %d;\n", idx, double(getView(view, "VectorType", 4)));
    pixels = double(getView(view, "ArrowSizePixels", 0));
    if pixels > 0
        fprintf(fid, "View[%d].ArrowSizeMin = %.17g;\n", idx, pixels);
        fprintf(fid, "View[%d].ArrowSizeMax = %.17g;\n", idx, pixels);
    end
elseif kind == "displacement"
    fprintf(fid, "View[%d].VectorType = 5;\n", idx);
    fprintf(fid, "View[%d].DisplacementFactor = %.17g;\n", ...
        idx, double(getView(view, "DisplacementFactor", 1.0)));
    fprintf(fid, "View[%d].DrawSkinOnly = %d;\n", idx, double(getView(view, "DrawSkinOnly", true)));
    fprintf(fid, "View[%d].SmoothNormals = %d;\n", idx, double(getView(view, "SmoothNormals", true)));
end
end


function value = getView(view, name, defaultValue)
if isfield(view, name)
    value = view.(name);
else
    value = defaultValue;
end
end


function rotation = cameraRotation(preset, customRotation)
if all(isfinite(customRotation))
    rotation = customRotation;
    return
end
switch preset
    case "z_up_xz_from_positive_y"
        rotation = [-68 0 0];
    case "positive_y_oblique"
        rotation = [75 0 -20];
    case "front_xz"
        rotation = [-90 0 0];
    otherwise
        rotation = [0 0 0];
end
end


function target = gmshMergeTarget(mshPath, geoPath)
[mshFolder, mshBase, mshExt] = fileparts(mshPath);
[geoFolder, ~, ~] = fileparts(geoPath);
if strcmpi(char(mshFolder), char(geoFolder))
    target = mshBase + mshExt;
else
    target = mshPath;
end
target = string(strrep(char(target), '\', '/'));
end


function names = viewNames(views)
names = strings(numel(views), 1);
for k = 1:numel(views)
    names(k) = string(getView(views(k), "Name", "view_" + string(k - 1)));
end
end


function maxIndex = maxViewIndexFrom(views)
maxIndex = -1;
for k = 1:numel(views)
    maxIndex = max(maxIndex, double(getView(views(k), "Index", k - 1)));
end
end


function s = escapeGmshString(s)
s = string(strrep(char(s), '"', '\"'));
end


function writeString(path, text)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeGmshPostLaunchArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", char(text));
clear cleanup
end
