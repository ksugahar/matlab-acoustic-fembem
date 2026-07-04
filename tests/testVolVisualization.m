function tests = testVolVisualization
%TESTVOLVISUALIZATION Lightweight .vol preview and summary tests.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testVolMeshSummaryFixture(testCase)
report = acoustic_fembem.vol_mesh_summary("unit_tetra.vol");

verifyEqual(testCase, report.status, "ok");
verifyEqual(testCase, report.tool, "acoustic_fembem_vol_mesh_summary");
verifyEqual(testCase, report.points, 4);
verifyEqual(testCase, report.triangles, 4);
verifyEqual(testCase, report.tets, 1);
verifyEqual(testCase, report.recommended_gui_viewer, "Netgen/native .vol viewer");
verifyEqual(testCase, report.recommended_matlab_preview, "plotVolMesh");
end


function testPlotVolMeshReturnsPatch(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol");
fig = figure("Visible", "off");
testCase.addTeardown(@() close(fig));

h = plotVolMesh(volFile);

verifyClass(testCase, h, "matlab.graphics.primitive.Patch");
verifyEqual(testCase, size(h.Faces, 1), 4);
end
