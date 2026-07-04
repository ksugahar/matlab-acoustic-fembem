function tests = testVolFemBemCqGmshArtifact
%TESTVOLFEMBEMCQGMSHARTIFACT .vol P1 FEM/BEM CQ Gmsh artifact writer.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testWriterProducesP1FemBemGmshArtifact(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", ...
    "four_tet_interior_node.vol"));
outBase = string(fullfile("C:\temp", ...
    "test_vol_fembem_cq_gmsh_" + char(java.util.UUID.randomUUID())));

artifact = writeVolFemBemCqGmsh3dArtifact(volFile, ...
    "OutputBase", outBase, ...
    "NumTime", 8, ...
    "TimeStep", 0.04, ...
    "CouplingForm", "JohnsonNedelec", ...
    "QuadratureOrder", 1);

verifyEqual(testCase, artifact.schema, ...
    "cae-ai-lab.vol-fembem-cq-gmsh3d-build.v1");
verifyEqual(testCase, artifact.gmsh_msh_version, "4.1");
verifyEqual(testCase, artifact.coupling_form, "JohnsonNedelec");
verifyTrue(testCase, artifact.checks.vol_mesh_tri_tet);
verifyTrue(testCase, artifact.checks.p1_volume_fem);
verifyTrue(testCase, artifact.checks.p1_boundary_bem);
verifyTrue(testCase, artifact.checks.johnson_nedelec_calderon_form);
verifyTrue(testCase, artifact.checks.double_layer_k_included);
verifyTrue(testCase, artifact.checks.high_order_impedance_boundary_lane_recorded);
verifyTrue(testCase, artifact.checks.not_kelvin_boundary);
verifyTrue(testCase, artifact.checks.node_data_pressure_written);
verifyTrue(testCase, artifact.checks.node_data_boundary_density_written);
verifyLessThan(testCase, artifact.max_relative_residual, 1e-8);

mshPath = outBase + ".msh";
geoPath = outBase + ".geo";
optPath = outBase + ".opt";
geoOptPath = geoPath + ".opt";
mshOptPath = mshPath + ".opt";
matPath = outBase + ".mat";
jsonPath = outBase + ".build.json";
verifyTrue(testCase, isfile(mshPath));
verifyTrue(testCase, isfile(geoPath));
verifyTrue(testCase, isfile(optPath));
verifyTrue(testCase, isfile(geoOptPath));
verifyTrue(testCase, isfile(mshOptPath));
verifyTrue(testCase, isfile(matPath));
verifyTrue(testCase, isfile(jsonPath));
verifyEqual(testCase, artifact.gmsh_opt, optPath);
verifyEqual(testCase, artifact.gmsh_geo_opt, geoOptPath);
verifyEqual(testCase, artifact.gmsh_msh_opt, mshOptPath);
verifyTrue(testCase, artifact.checks.gmsh_geo_opt_exact_autoload_written);
verifyTrue(testCase, artifact.checks.gmsh_msh_opt_exact_autoload_written);

mshText = string(fileread(mshPath));
verifyTrue(testCase, contains(mshText, "$MeshFormat" + newline + "4.1 0 8"));
verifyTrue(testCase, contains(mshText, """boundary_p1_bem_surface"""));
verifyTrue(testCase, contains(mshText, """vol_p1_fem_acoustic_domain"""));
verifyTrue(testCase, contains(mshText, "$Elements" + newline + "2 "));
verifyTrue(testCase, contains(mshText, """interior_pressure"""));
verifyTrue(testCase, contains(mshText, """boundary_density"""));
verifyEqual(testCase, count(mshText, "$NodeData"), 16);

optText = string(fileread(optPath));
geoOptText = string(fileread(geoOptPath));
mshOptText = string(fileread(mshOptPath));
verifyTrue(testCase, contains(optText, "General.Trackball = 0"));
verifyTrue(testCase, contains(optText, "General.RotationX = 75"));
verifyTrue(testCase, contains(optText, "General.RotationY = 0"));
verifyTrue(testCase, contains(optText, "General.RotationZ = -20"));
verifyTrue(testCase, contains(optText, "Mesh.SurfaceEdges = 1"));
verifyTrue(testCase, contains(optText, "View[0].ColormapNumber = 2"));
verifyTrue(testCase, contains(optText, "Mesh.Clip = 0"));
verifyTrue(testCase, contains(optText, "View[0].Clip = 1"));
verifyTrue(testCase, contains(optText, "View[0].IntervalsType = 2"));
verifyEqual(testCase, geoOptText, optText);
verifyTrue(testCase, contains(mshOptText, "Gmsh mesh-inspection options"));
verifyTrue(testCase, contains(mshOptText, "General.InitialModule = 1"));
verifyTrue(testCase, contains(mshOptText, "Mesh.SurfaceFaces = 1"));
verifyTrue(testCase, contains(mshOptText, "View[0].Visible = 0"));

jsonText = string(fileread(jsonPath));
verifyTrue(testCase, startsWith(jsonText, "{""schema"":"));
verifyTrue(testCase, contains(jsonText, """p1_volume_fem"":true"));
verifyTrue(testCase, contains(jsonText, """p1_boundary_bem"":true"));
end
