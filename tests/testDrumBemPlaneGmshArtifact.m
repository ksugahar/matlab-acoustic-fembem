function tests = testDrumBemPlaneGmshArtifact
%TESTDRUMBEMPLANEGMSHARTIFACT Gmsh post artifact for BEM drum scene.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testWriterProducesPlanePressureAndDeformingDrumSurface(testCase)
scene = drumFemBemCoupledDemo( ...
    "NumX", 32, ...
    "NumZ", 32, ...
    "NumTime", 12, ...
    "TMax", 1.2e-3);
scene.status = "ok";  % this test targets the writer contract, not the solver gate.
src = string(fullfile("C:\temp", ...
    "test_drum_bem_scene_" + char(java.util.UUID.randomUUID()) + ".mat"));
save(src, "scene");
outBase = string(fullfile("C:\temp", ...
    "test_drum_bem_gmsh_" + char(java.util.UUID.randomUUID())));

artifact = writeDrumBemPlaneGmshArtifact(src, ...
    "OutputBase", outBase, ...
    "PlaneStride", 2, ...
    "DrumRadial", 4, ...
    "DrumAzimuth", 12, ...
    "DrumAxial", 4);

verifyEqual(testCase, artifact.schema, "cae-ai-lab.drum-bem-plane-gmsh.v1");
verifyTrue(testCase, artifact.pass);
verifyEqual(testCase, artifact.acoustic_view, "x_z_plane_scalar_pressure");
verifyEqual(testCase, artifact.drum_view, "top_bottom_side_bem_surfaces_vector_displacement");
verifyTrue(testCase, artifact.checks.pressure_xz_plane_written);
verifyTrue(testCase, artifact.checks.drum_3d_surface_written);
verifyTrue(testCase, artifact.checks.deforming_surface_written);
verifyGreaterThan(testCase, artifact.gmsh_displacement_factor, 1);

mshPath = outBase + ".msh";
geoPath = outBase + ".geo";
geoOptPath = geoPath + ".opt";
mshOptPath = mshPath + ".opt";
jsonPath = outBase + ".result.json";
verifyTrue(testCase, isfile(mshPath));
verifyTrue(testCase, isfile(geoPath));
verifyTrue(testCase, isfile(geoOptPath));
verifyTrue(testCase, isfile(mshOptPath));
verifyTrue(testCase, isfile(jsonPath));

mshText = string(fileread(mshPath));
verifyTrue(testCase, contains(mshText, """bem_pressure_xz_plane"""));
verifyTrue(testCase, contains(mshText, """drum_top_head_bem_surface"""));
verifyTrue(testCase, contains(mshText, """drum_bottom_head_bem_surface"""));
verifyTrue(testCase, contains(mshText, """drum_side_shell_bem_surface"""));
verifyTrue(testCase, contains(mshText, """drum_surface_displacement"""));
verifyEqual(testCase, count(mshText, "$NodeData"), 24);

geoText = string(fileread(geoPath));
geoOptText = string(fileread(geoOptPath));
verifyTrue(testCase, contains(geoText, "BEM acoustic field on x-z plane"));
verifyTrue(testCase, contains(geoOptText, "View[0].Name = ""BEM pressure on x-z plane"""));
verifyTrue(testCase, contains(geoOptText, "View[1].VectorType = 5"));
verifyTrue(testCase, contains(geoOptText, "View[1].DisplacementFactor"));
verifyTrue(testCase, contains(geoOptText, "PostProcessing.Link = 1"));
verifyTrue(testCase, contains(geoOptText, "PostProcessing.AnimationCycle = 0"));
verifyTrue(testCase, contains(geoOptText, "General.RotationX = -68"));
verifyTrue(testCase, contains(geoOptText, "General.RotationY = 0"));
verifyTrue(testCase, contains(geoOptText, "General.RotationZ = 0"));

mshOptText = string(fileread(mshOptPath));
verifyTrue(testCase, contains(mshOptText, "Open .geo for the post display"));
verifyTrue(testCase, contains(mshOptText, "Mesh.SurfaceFaces = 1"));
verifyTrue(testCase, contains(mshOptText, "View[0].Visible = 0"));
verifyTrue(testCase, contains(mshOptText, "General.RotationX = -68"));
end
