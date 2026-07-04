function tests = testDrumStepTimeField
%TESTDRUMSTEPTIMEFIELD Time-domain drum/Rayleigh teaching example.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testStepFieldIsCausalAndFinite(testCase)
field = drumStepTimeField( ...
    "NumRadialObservation", 18, ...
    "NumAxialObservation", 16, ...
    "NumSourceRadial", 8, ...
    "NumSourceAzimuth", 16, ...
    "NumTime", 24, ...
    "TMax", 8e-4);

verifyEqual(testCase, field.kind, "drum_step_time_field_rayleigh");
verifySize(testCase, field.pressure, [18, 16, 24]);
verifyTrue(testCase, field.checks.causal_initial_field_zero);
verifyTrue(testCase, field.checks.finite_pressure);
verifyTrue(testCase, field.checks.nonzero_after_wave_arrival);
verifyGreaterThan(testCase, field.summary.max_abs_pressure, 0);
end


function testStepFieldPlotReturnsAxes(testCase)
field = drumStepTimeField( ...
    "NumRadialObservation", 12, ...
    "NumAxialObservation", 10, ...
    "NumSourceRadial", 6, ...
    "NumSourceAzimuth", 12, ...
    "NumTime", 12, ...
    "TMax", 6e-4);

fig = figure("Visible", "off");
testCase.addTeardown(@() close(fig));
ax = axes(fig);
out = plotDrumStepTimeField(field, 6, "Parent", ax);

verifyEqual(testCase, out, ax);
verifyEqual(testCase, numel(ax.Children), 1);
end


function testStepFieldGifWritesAnimation(testCase)
field = drumStepTimeField( ...
    "NumRadialObservation", 14, ...
    "NumAxialObservation", 12, ...
    "NumSourceRadial", 6, ...
    "NumSourceAzimuth", 12, ...
    "NumTime", 10, ...
    "TMax", 6e-4);

outDir = "C:\temp";
if ~isfolder(outDir)
    mkdir(outDir);
end
gifPath = string(tempname(outDir)) + ".gif";
testCase.addTeardown(@() deleteIfExists(gifPath));

info = writeDrumStepTimeGif(field, gifPath, ...
    "TimeIndices", [1, 4, 7, 10], ...
    "DelayTime", 0.01);

verifyEqual(testCase, info.kind, "drum_step_time_field_gif");
verifyEqual(testCase, info.num_frames, 4);
verifyTrue(testCase, isfile(gifPath));

frames = imfinfo(gifPath);
verifyEqual(testCase, numel(frames), 4);
verifyGreaterThan(testCase, frames(1).Width, 1);
verifyGreaterThan(testCase, frames(1).Height, 1);
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
