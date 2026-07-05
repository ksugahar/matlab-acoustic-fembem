function tests = testGmshPostLaunchArtifact
%TESTGMSHPOSTLAUNCHARTIFACT Generic Gmsh .geo/.opt launch writer.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testWriterRecordsCameraCutPlaneAndViews(testCase)
outBase = string(fullfile("C:\temp", ...
    "test_gmsh_post_launch_" + char(java.util.UUID.randomUUID())));
mshPath = outBase + ".msh";
writeText(mshPath, "$MeshFormat" + newline + "4.1 0 8" + newline + "$EndMeshFormat");

views = defaultViewStruct();
views(1).Index = 0;
views(1).Name = "BEM pressure on x-z plane";
views(1).Kind = "scalar";
views(1).TimeStep = 10;
views(1).Range = [-2.0 2.0];
views(1).NbIso = 40;
views(1).ColormapNumber = 2;
views(1).ShowScale = true;

views(2) = defaultViewStruct();
views(2).Index = 1;
views(2).Name = "3D deforming drum BEM surface";
views(2).Kind = "displacement";
views(2).TimeStep = 10;
views(2).ShowElement = true;
views(2).ShowScale = false;
views(2).ColormapNumber = 10;
views(2).DisplacementFactor = 3.5;

artifact = writeGmshPostLaunchArtifact(mshPath, ...
    "OutputBase", outBase, ...
    "Title", "Unit-test Gmsh post launch", ...
    "CameraPreset", "z_up_xz_from_positive_y", ...
    "EnableCutPlane", true, ...
    "CutPlaneNormal", [0 -1 0], ...
    "CutPlaneOffset", 0, ...
    "MeshSurfaceFaces", false, ...
    "MeshSurfaceEdges", true, ...
    "Views", views);

verifyEqual(testCase, artifact.schema, "cae-ai-lab.gmsh-post-launch.v1");
verifyTrue(testCase, artifact.pass);
verifyEqual(testCase, artifact.launch_target, outBase + ".geo");
verifyEqual(testCase, artifact.mesh_inspection_target, mshPath);
verifyEqual(testCase, artifact.camera_preset, "z_up_xz_from_positive_y");
verifyEqual(testCase, artifact.cut_plane_enabled, true);
verifyEqual(testCase, artifact.view_count, 2);
verifyEqual(testCase, artifact.view_names(1), "BEM pressure on x-z plane");

geoPath = outBase + ".geo";
geoOptPath = geoPath + ".opt";
optPath = outBase + ".opt";
mshOptPath = mshPath + ".opt";
jsonPath = outBase + ".display.json";
verifyTrue(testCase, isfile(geoPath));
verifyTrue(testCase, isfile(geoOptPath));
verifyTrue(testCase, isfile(optPath));
verifyTrue(testCase, isfile(mshOptPath));
verifyTrue(testCase, isfile(jsonPath));

geoText = string(fileread(geoPath));
geoOptText = string(fileread(geoOptPath));
mshOptText = string(fileread(mshOptPath));
verifyTrue(testCase, contains(geoText, "Merge ""test_gmsh_post_launch_"));
verifyTrue(testCase, contains(geoText, "General.RotationX = -68"));
verifyTrue(testCase, contains(geoOptText, "General.Clip0B = -1"));
verifyTrue(testCase, contains(geoOptText, "Mesh.Clip = 1"));
verifyTrue(testCase, contains(geoOptText, "View[0].Name = ""BEM pressure on x-z plane"""));
verifyTrue(testCase, contains(geoOptText, "View[0].CustomMin = -2"));
verifyTrue(testCase, contains(geoOptText, "View[1].VectorType = 5"));
verifyTrue(testCase, contains(geoOptText, "View[1].DisplacementFactor = 3.5"));
verifyTrue(testCase, contains(geoOptText, "PostProcessing.AnimationCycle = 0"));
verifyTrue(testCase, contains(mshOptText, "Open .geo for the post display"));
verifyTrue(testCase, contains(mshOptText, "View[0].Visible = 0"));
verifyTrue(testCase, contains(mshOptText, "View[1].Visible = 0"));

jsonText = string(fileread(jsonPath));
verifyTrue(testCase, contains(jsonText, """schema"":""cae-ai-lab.gmsh-post-launch.v1"""));
verifyTrue(testCase, contains(jsonText, """launch_target_is_geo"":true"));
end


function v = defaultViewStruct()
v = struct( ...
    "Index", 0, ...
    "Name", "", ...
    "Kind", "scalar", ...
    "Visible", true, ...
    "TimeStep", 0, ...
    "IntervalsType", 2, ...
    "NbIso", 32, ...
    "ColormapNumber", 2, ...
    "ShowElement", false, ...
    "ShowScale", true, ...
    "Light", false, ...
    "Clip", true, ...
    "Range", [], ...
    "VectorType", 4, ...
    "ArrowSizePixels", 0, ...
    "DisplacementFactor", 1.0, ...
    "DrawSkinOnly", true, ...
    "SmoothNormals", true);
end


function writeText(path, text)
fid = fopen(path, "w", "n", "UTF-8");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", char(text));
clear cleanup
end
