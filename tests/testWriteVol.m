function tests = testWriteVol
%TESTWRITEVOL Pure-MATLAB .vol writing: writeVol + icosphereVol round-trips.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testWriteVolTetRoundTripsThroughVolMesh(testCase)
path = tempPath(testCase);
points = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
tets = [1 2 3 4];
tri = [1 3 2; 1 2 4; 2 3 4; 1 4 3];

report = writeVol(path, points, tri, Tets=tets);
verifyEqual(testCase, report.status, "ok");
verifyFalse(testCase, report.surface_only);

mesh = VolMesh(path);
verifyEqual(testCase, mesh.vtx, points, "AbsTol", 1e-12);
verifyEqual(testCase, size(mesh.tri, 1), 4);
verifyEqual(testCase, size(mesh.tet, 1), 1);
verifyEqual(testCase, string(mesh.materials(1)), "domain");
verifyEqual(testCase, string(mesh.boundaries(1)), "outer");
end


function testWriteVolSurfaceOnlyIsBemReady(testCase)
% no tets: the boundary triangles alone must still drive a SurfaceMesh
path = tempPath(testCase);
points = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
tri = [1 3 2; 1 2 4; 2 3 4; 1 4 3];

report = writeVol(path, points, tri);           % Tets omitted
verifyTrue(testCase, report.surface_only);
verifyEqual(testCase, report.tets, 0);

mesh = VolMesh(path);
verifyEqual(testCase, size(mesh.tet, 1), 0);
surface = mesh.boundary();
verifyEqual(testCase, size(surface.tri, 1), 4);
verifyEqual(testCase, numel(surface.areas()), 4);
end


function testWriteVolNamesMultipleBoundaries(testCase)
path = tempPath(testCase);
points = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
tri = [1 3 2; 1 2 4; 2 3 4; 1 4 3];

writeVol(path, points, tri, ...
    TriBoundaryId=[1; 1; 2; 2], BoundaryNames=["floor" "roof"]);
mesh = VolMesh(path);
verifyEqual(testCase, string(mesh.boundaries(1)), "floor");
verifyEqual(testCase, string(mesh.boundaries(2)), "roof");
end


function testIcosphereVolIsClosedOutwardUnitSphere(testCase)
path = tempPath(testCase);
report = icosphereVol(path, Subdivisions=3);

verifyEqual(testCase, report.status, "ok");
verifyEqual(testCase, report.surface_nodes, 642);       % 10*4^3 + 2
verifyEqual(testCase, report.boundary_orientation, "outward");  % cone tets present

mesh = VolMesh(path);
surface = mesh.boundary();
verifyEqual(testCase, size(surface.vtx, 1), 642);
verifyEqual(testCase, size(surface.tri, 1), 1280);      % 20*4^3

% flat triangles chord the sphere: area just under 4*pi, well within 2%
area = sum(surface.areas());
verifyLessThan(testCase, area, 4*pi + 1e-9);
verifyGreaterThan(testCase, area, 4*pi * 0.98);

% every vertex sits on the unit sphere
radii = vecnorm(surface.vtx, 2, 2);
verifyEqual(testCase, radii, ones(size(radii)), "AbsTol", 1e-12);
end


function testIcosphereVolSurfaceOnlyOption(testCase)
path = tempPath(testCase);
report = icosphereVol(path, Subdivisions=2, SurfaceOnly=true);
verifyTrue(testCase, report.surface_only);

mesh = VolMesh(path);
verifyEqual(testCase, size(mesh.tet, 1), 0);
verifyEqual(testCase, size(mesh.boundary().tri, 1), 320);   % 20*4^2
end


function path = tempPath(testCase)
path = string(fullfile(tempdir, "writeVolTest_" + ...
    char(java.util.UUID.randomUUID()) + ".vol"));
testCase.addTeardown(@() deleteIfExists(path));
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
