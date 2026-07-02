function tests = testReadVolTriTet
%TESTREADVOLTRITET Unit tests for the MATLAB .vol tri/tet intake.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testReadTetVol(testCase)
path = writeFixture(testCase, tetVolText());

mesh = readVolTriTet(path);

verifyEqual(testCase, size(mesh.vtx), [4 3]);
verifyEqual(testCase, mesh.tri, [1 2 3; 1 4 2; 2 4 3; 3 4 1]);
verifyEqual(testCase, mesh.triCol, ones(4, 1));
verifyEqual(testCase, mesh.tet, [1 2 3 4]);
verifyEqual(testCase, mesh.tetMat, 1);
verifyEqual(testCase, mesh.materials(1), 'air');
verifyEqual(testCase, mesh.boundaries(1), 'outer');
verifyEqual(testCase, mesh.traceNodeIds, (1:4).');
verifyEqual(testCase, mesh.policy, "netgen_vol_tri_tet_only_shared_one_based_nodes");
end


function testBoundaryOrientationSummaryUsesAdjacentTet(testCase)
path = writeFixture(testCase, tetVolText());

mesh = readVolTriTet(path);

verifyEqual(testCase, mesh.boundaryOrientation.boundaryOrientation, "inward");
verifyEqual(testCase, mesh.boundaryOrientation.triangleOrientationSignsToOutward, -ones(4, 1));
verifyEqual(testCase, mesh.boundaryOrientation.adjacentTetIndices, ones(4, 1));
verifyEqual(testCase, mesh.boundaryOrientation.rows(1).normalOrientation, "inward");
verifyEqual(testCase, mesh.boundaryOrientation.rows(1).storedAreaVector, [0 0 0.5], "AbsTol", 1e-12);
verifyEqual(testCase, mesh.boundaryOrientation.rows(1).outwardAreaVector, [0 0 -0.5], "AbsTol", 1e-12);
end


function testReadSurfaceElementsUvTriVol(testCase)
path = writeFixture(testCase, surfaceElementsUvTetVolText());

mesh = readVolTriTet(path);

verifyEqual(testCase, size(mesh.vtx), [4 3]);
verifyEqual(testCase, mesh.tri, [1 2 3; 1 4 2; 2 4 3; 3 4 1]);
verifyEqual(testCase, mesh.tet, [1 2 3 4]);
verifyEqual(testCase, mesh.summary.triangles, 4);
verifyEqual(testCase, mesh.summary.tets, 1);
end


function testFemBemModelUsesSharedOneBasedNodes(testCase)
path = writeFixture(testCase, tetVolText());

model = FemBemModel(path);

verifyEqual(testCase, model.mesh.tet, [1 2 3 4]);
verifyEqual(testCase, model.surface.tri(1, :), [1 2 3]);
verifyEqual(testCase, model.surface.col, ones(4, 1));
verifyEqual(testCase, model.mesh.traceNodeIds, (1:4).');
verifyEqual(testCase, model.surface.orientation.boundaryOrientation, "inward");
verifyEqual(testCase, model.surface.orientation.triangleOrientationSignsToOutward, -ones(4, 1));
verifyEqual(testCase, model.surface.orientation.adjacentTetIndices, ones(4, 1));
verifyEqual(testCase, [model.trace.rowIdentity.trace_row_index].', (1:4).');
verifyEqual(testCase, [model.trace.rowIdentity.fem_node_id].', (1:4).');
verifyEqual(testCase, [model.trace.rowIdentity.bem_node_id].', (1:4).');
verifyEqual(testCase, model.status, "vol_ready_first_order_h1_hcurl_rwg");
end


function testRejectQuadSurface(testCase)
bad = replace(tetVolText(), "1 1 1 0 3 1 2 3", "1 1 1 0 4 1 2 3 4");
path = writeFixture(testCase, bad);

verifyError(testCase, @() readVolTriTet(path), "readVolTriTet:surface");
end


function testRejectQuadSurfaceElementsUv(testCase)
bad = replace(surfaceElementsUvTetVolText(), ...
    "1 1 1 0 3 1 2 3 0 0 1 0 0 1", ...
    "1 1 1 0 4 1 2 3 4 0 0 1 0 1 1 0 1");
path = writeFixture(testCase, bad);

verifyError(testCase, @() readVolTriTet(path), "readVolTriTet:surface");
end


function testRejectHexVolume(testCase)
bad = replace(tetVolText(), "1 4 1 2 3 4", "1 8 1 2 3 4 1 2 3 4");
path = writeFixture(testCase, bad);

verifyError(testCase, @() readVolTriTet(path), "readVolTriTet:volume");
end


function testRejectCurvedElements(testCase)
bad = replace(tetVolText(), "endmesh", join([
    "curvedelements"
    "1"
    "1 2 3 4 5"
    "endmesh"
    ], newline));
path = writeFixture(testCase, bad);

verifyError(testCase, @() readVolTriTet(path), "readVolTriTet:curved");
end


function path = writeFixture(testCase, text)
path = string(fullfile(tempdir, "readVolTriTet_fixture_" + char(java.util.UUID.randomUUID()) + ".vol"));
fid = fopen(path, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", text);
clear cleanup
testCase.addTeardown(@() delete(path));
end


function text = tetVolText()
text = join([
    "mesh3d"
    "dimension"
    "3"
    "geomtype"
    "0"
    "facedescriptors"
    "1"
    "1 1 0 1 1"
    "surfaceelements"
    "4"
    "1 1 1 0 3 1 2 3"
    "1 1 1 0 3 1 4 2"
    "1 1 1 0 3 2 4 3"
    "1 1 1 0 3 3 4 1"
    "volumeelements"
    "1"
    "1 4 1 2 3 4"
    "points"
    "4"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "pointelements"
    "0"
    "materials"
    "1"
    "1 air"
    "bcnames"
    "1"
    "1 outer"
    "endmesh"
    ], newline);
end


function text = surfaceElementsUvTetVolText()
text = join([
    "mesh3d"
    "dimension"
    "3"
    "geomtype"
    "0"
    "facedescriptors"
    "1"
    "1 1 0 1 1"
    "surfaceelementsuv"
    "4"
    "1 1 1 0 3 1 2 3 0 0 1 0 0 1"
    "1 1 1 0 3 1 4 2 0 0 0 1 1 0"
    "1 1 1 0 3 2 4 3 1 0 0 1 0 0"
    "1 1 1 0 3 3 4 1 0 1 1 0 0 0"
    "volumeelements"
    "1"
    "1 4 1 2 3 4"
    "points"
    "4"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "pointelements"
    "0"
    "materials"
    "1"
    "1 air"
    "bcnames"
    "1"
    "1 outer"
    "endmesh"
    ], newline);
end
