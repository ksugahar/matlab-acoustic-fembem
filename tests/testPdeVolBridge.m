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


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
