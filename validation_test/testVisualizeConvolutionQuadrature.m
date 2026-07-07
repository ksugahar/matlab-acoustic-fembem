function tests = testVisualizeConvolutionQuadrature
%TESTVISUALIZECONVOLUTIONQUADRATURE Smoke tests for the CQ solver X-ray figure.
%
%   visualizeConvolutionQuadrature turns the internals already returned by
%   volTdBemConvolutionQuadrature into a six-panel figure.  These tests confirm
%   it builds the figure (from a result struct and from a .vol path) and honours
%   the SavePath export, without asserting pixels.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testBuildsSixPanelFigureFromResult(testCase)
vol = sphereFixture();
result = volTdBemConvolutionQuadrature(vol, NumTime=16, TimeStep=0.3, ...
    Method="BDF2", QuadratureOrder=1);
png = tempPng();
testCase.addTeardown(@() deleteIfExists(png));

fig = visualizeConvolutionQuadrature(result, SavePath=png, Verbose=false);
testCase.addTeardown(@() close(fig));

verifyTrue(testCase, isgraphics(fig, "figure"));
ax = findobj(fig, "Type", "axes");
verifyGreaterThanOrEqual(testCase, numel(ax), 6);   % six method panels
verifyTrue(testCase, isfile(png));                  % SavePath exported a PNG
end


function testRunsSolverFromVolPath(testCase)
% Passing a .vol (not a struct) runs the solver internally first.
vol = sphereFixture();
png = tempPng();
testCase.addTeardown(@() deleteIfExists(png));

fig = visualizeConvolutionQuadrature(vol, NumTime=12, TimeStep=0.3, ...
    Method="BDF1", QuadratureOrder=1, SavePath=png, Verbose=false);
testCase.addTeardown(@() close(fig));

verifyTrue(testCase, isgraphics(fig, "figure"));
verifyTrue(testCase, isfile(png));
end


function vol = sphereFixture()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
vol = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_coarse.vol"));
end


function p = tempPng()
p = string(fullfile(tempdir, "cq_xray_" + char(java.util.UUID.randomUUID()) + ".png"));
end


function deleteIfExists(p)
if isfile(p)
    delete(p);
end
end
