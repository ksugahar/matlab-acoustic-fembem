function tests = testReadVolTriTet
%TESTREADVOLTRITET Unit tests for the MATLAB .vol tri/tet intake.

tests = functiontests(localfunctions);
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


function testVolFemBemModelUsesSharedOneBasedNodes(testCase)
path = writeFixture(testCase, tetVolText());

model = volFemBemModel(path);

verifyEqual(testCase, model.lukas.geo.conn_matrix, [1 2 3 4]);
verifyEqual(testCase, model.gypsilab.elt(1, :), [1 2 3]);
verifyEqual(testCase, model.gypsilab.col, ones(4, 1));
verifyEqual(testCase, model.trace.nodeIds, (1:4).');
verifyEqual(testCase, model.trace.boundaryOrientation, "inward");
verifyEqual(testCase, model.trace.triangleOrientationSignsToOutward, -ones(4, 1));
verifyEqual(testCase, model.trace.adjacentTetIndices, ones(4, 1));
verifyEqual(testCase, [model.trace.traceRowIdentity.trace_row_index].', (1:4).');
verifyEqual(testCase, [model.trace.traceRowIdentity.fem_node_id].', (1:4).');
verifyEqual(testCase, [model.trace.traceRowIdentity.bem_node_id].', (1:4).');
verifyEqual(testCase, model.status, "mesh_ready");
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
