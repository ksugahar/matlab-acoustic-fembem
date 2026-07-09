function tests = testCurvedPanelNetgenSmoke
%TESTCURVEDPANELNETGENSMOKE Fast smoke for the netgen curved-node companion path.
%
% Reads the committed netgen fixture (fixtures/curved_panels/ng_sphere.vol +
% ng_sphere_curved_p2.json), matches the curved nodes onto the SurfaceMesh, and
% checks the curved-panel area beats the flat one -- the heavy accuracy A/B is
% validation_test/testCurvedPanelNetgenVol.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testNetgenCompanionBuildsCurvedQuadrature(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
vol = fullfile(repoRoot, "fixtures", "curved_panels", "ng_sphere.vol");
json = fullfile(repoRoot, "fixtures", "curved_panels", "ng_sphere_curved_p2.json");
surface = VolMesh(vol).boundary();
R = mean(vecnorm(surface.vtx, 2, 2));

qN = curvedQuadratureFromNetgen(surface, json, 7);
verifyEqual(testCase, qN.curveOrder, 2);
verifyEqual(testCase, size(qN.geomNodes), [size(surface.tri,1), 6, 3]);
verifyTrue(testCase, all(isfinite(qN.geomNodes(:))));

areaTrue = 4 * pi * R^2;
areaFlat = sum(CurvedPanelQuadrature(surface, 7, @(X) X, 1).weights);
areaN = sum(qN.weights);
verifyLessThan(testCase, abs(areaN - areaTrue) / areaTrue, 5e-3);
verifyLessThan(testCase, abs(areaN - areaTrue), abs(areaFlat - areaTrue));
end
