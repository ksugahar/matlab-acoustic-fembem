function tests = testPdeVolBridge
%TESTPDEVOLBRIDGE PDE Toolbox mesh-to-.vol bridge tests.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testSyntheticPdeMeshWritesReadableVol(testCase)
pdeMesh = struct();
pdeMesh.Nodes = [
    0 1 0 0
    0 0 1 0
    0 0 0 1
];
pdeMesh.Elements = [1; 2; 3; 4];

volFile = string(fullfile(tempdir, "pde_mesh_bridge_" + char(java.util.UUID.randomUUID()) + ".vol"));
testCase.addTeardown(@() deleteIfExists(volFile));

report = writePdeMeshVol(pdeMesh, volFile, ...
    MaterialName="air", BoundaryName="outer");
mesh = readVolTriTet(volFile);

verifyEqual(testCase, report.status, "ok");
verifyEqual(testCase, report.policy, "pde_toolbox_linear_tet_to_netgen_vol_tri_tet");
verifyEqual(testCase, report.points, 4);
verifyEqual(testCase, report.triangles, 4);
verifyEqual(testCase, report.tets, 1);
verifyEqual(testCase, mesh.summary.points, 4);
verifyEqual(testCase, mesh.summary.triangles, 4);
verifyEqual(testCase, mesh.summary.tets, 1);
verifyEqual(testCase, mesh.materials(1), 'air');
verifyEqual(testCase, mesh.boundaries(1), 'outer');
verifyEqual(testCase, mesh.boundaryOrientation.boundaryOrientation, "outward");
end


function testRejectsNonLinearTetMesh(testCase)
pdeMesh = struct();
pdeMesh.Nodes = zeros(3, 10);
pdeMesh.Elements = (1:10).';

verifyError(testCase, ...
    @() writePdeMeshVol(pdeMesh, fullfile(tempdir, "bad_pde_mesh.vol")), ...
    "writePdeMeshVol:elementOrder");
end


function testWritePdeBoxVol(testCase)
% writePdeBoxVol needs the OPTIONAL PDE Toolbox.  WITH it, it meshes a box;
% WITHOUT it, it must fail loud - both branches PASS, so run_tests stays green on
% a minimal (no-toolbox) checkout too.  Also guards the R2026a built-in-multicuboid
% regression (the availability check must accept a built-in, not only a .m file).
volFile = string(fullfile(tempdir, "pde_box_" + char(java.util.UUID.randomUUID()) + ".vol"));
testCase.addTeardown(@() deleteIfExists(volFile));

if exist("createpde") == 0 || exist("multicuboid") == 0
    verifyError(testCase, ...
        @() writePdeBoxVol(volFile, Size=[1.2 1.0 0.8], Hmax=0.4), ...
        "writePdeBoxVol:pdeToolboxUnavailable");
    return
end

report = writePdeBoxVol(volFile, Size=[1.2 1.0 0.8], Hmax=0.4, ...
    MaterialName="solid", BoundaryName="interface");
verifyEqual(testCase, report.status, "ok");
verifyGreaterThan(testCase, report.tets, 0);
verifyGreaterThan(testCase, report.triangles, 0);
verifyEqual(testCase, report.boundary_orientation, "outward");

mesh = VolMesh(volFile);
verifyEqual(testCase, mesh.summary.tets, report.tets);
end


function testWritePdeGeometryVolSphere(testCase)
% The general geometry->.vol path meshes ANY PDE Toolbox geometry (here a sphere
% via multisphere), not just a box.  WITH the toolbox it meshes; WITHOUT it,
% writePdeGeometryVol fails loud before touching the geometry - both branches PASS.
volFile = string(fullfile(tempdir, "pde_sphere_" + char(java.util.UUID.randomUUID()) + ".vol"));
testCase.addTeardown(@() deleteIfExists(volFile));

if exist("createpde") == 0
    verifyError(testCase, ...
        @() writePdeGeometryVol(volFile, [], Hmax=0.15), ...
        "writePdeGeometryVol:pdeToolboxUnavailable");
    return
end

report = writePdeGeometryVol(volFile, multisphere(0.5), Hmax=0.15, ...
    MaterialName="solid", BoundaryName="interface");
verifyEqual(testCase, report.status, "ok");
verifyGreaterThan(testCase, report.tets, 0);
verifyEqual(testCase, report.boundary_orientation, "outward");
verifyEqual(testCase, report.generator, "pde_toolbox_geometry");
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
