function tests = testFirstOrderFemBemSpaces
%TESTFIRSTORDERFEMBEMSPACES Tests for Gypsilab-style .vol first-order spaces.

tests = functiontests(localfunctions);
end


function testSimpleApiBuildsH1HcurlRwg(testCase)
path = writeFixture(testCase, tetVolText());

m = volFemBem(path);

verifyEqual(testCase, m.status, "vol_ready_first_order_h1_hcurl_rwg");
verifyEqual(testCase, m.h1.basis, "P1");
verifyEqual(testCase, m.hcurl.basis, "Nedelec0");
verifyEqual(testCase, m.rwg.basis, "RWG0");
verifyEqual(testCase, size(m.hcurl.edges, 1), 6);
verifyEqual(testCase, numel(m.rwg.dofEdgeIds), 6);
verifyEqual(testCase, m.rwg.hcurlEdgeIds(:), (1:6).');
end


function testH1StiffnessOnUnitTetra(testCase)
path = writeFixture(testCase, tetVolText());
m = volFemBem(path);

fem = m.h1.stiffness();
K = full(fem.stiffness);
expected = (1 / 6) * [ ...
     3 -1 -1 -1
    -1  1  0  0
    -1  0  1  0
    -1  0  0  1];

verifyEqual(testCase, K, expected, "AbsTol", 1e-14);
verifyEqual(testCase, sum(K, 2), zeros(4, 1), "AbsTol", 1e-14);
verifyEqual(testCase, fem.volumes, 1/6, "AbsTol", 1e-14);
end


function testBoundaryCompactionKeepsTraceToInteriorVolumeNodes(testCase)
path = writeFixture(testCase, fourTetWithInteriorNodeVolText());

m = volFemBem(path);

verifyEqual(testCase, size(m.h1.nodes, 1), 5);
verifyEqual(testCase, size(m.rwg.vtx, 1), 4);
verifyEqual(testCase, m.rwg.globalNodeIds, (1:4).');
verifyEqual(testCase, size(m.topology.trace.h1ToScalarBem), [4 5]);

u = (10:10:50).';
verifyEqual(testCase, m.topology.trace.h1ToScalarBem * u, u(1:4));
end


function testAssembledTraceScaffoldContainsHcurlRwgMap(testCase)
path = writeFixture(testCase, tetVolText());
m = volFemBemModel(path);

m = assembleFirstOrderFemBemTrace(m);

verifyEqual(testCase, m.status, "operators_ready_first_order_h1_hcurl_rwg_trace");
verifyEqual(testCase, m.operators.trace.rwgToHcurlEdgeIds(:), (1:6).');
verifyEqual(testCase, size(m.operators.fem.stiffness), [4 4]);
verifyEqual(testCase, size(m.operators.bem.surfaceMass), [4 4]);
end


function path = writeFixture(testCase, text)
path = string(fullfile(tempdir, "firstOrderFemBem_" + char(java.util.UUID.randomUUID()) + ".vol"));
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


function text = fourTetWithInteriorNodeVolText()
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
    "1 1 1 0 3 1 3 4"
    "1 1 1 0 3 2 4 3"
    "volumeelements"
    "4"
    "1 4 1 2 3 5"
    "1 4 1 4 2 5"
    "1 4 1 3 4 5"
    "1 4 2 4 3 5"
    "points"
    "5"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "0.25 0.25 0.25"
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
