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


function testHighOrderImpedanceSceneAndGif(testCase)
field = drumStepTimeField( ...
    "NumRadialObservation", 24, ...
    "NumAxialObservation", 22, ...
    "NumSourceRadial", 7, ...
    "NumSourceAzimuth", 14, ...
    "NumTime", 12, ...
    "TMax", 7e-4);

scene = drumHighOrderImpedanceScene(field, ...
    "NumX", 50, ...
    "NumZ", 42, ...
    "TimeIndices", [1, 5, 9, 12]);

verifyEqual(testCase, scene.kind, "drum_high_order_impedance_scene");
verifyTrue(testCase, contains(scene.boundary_type, "high-order impedance"));
verifySize(testCase, scene.pressure, [50, 50, 4]);
verifyTrue(testCase, scene.axis.equal);
verifyEqual(testCase, scene.x(1), scene.z(1), "AbsTol", 1e-12);
verifyEqual(testCase, scene.x(end), scene.z(end), "AbsTol", 1e-12);
verifyEqual(testCase, scene.geometry.struck_surface, "top membrane at z=0");
verifyTrue(testCase, any(scene.masks.high_order_impedance_boundary, "all"));
verifyTrue(testCase, any(scene.masks.high_order_impedance_boundary(scene.z < 0, :), "all"));
verifyTrue(testCase, any(scene.masks.drum_frame, "all"));
verifyTrue(testCase, any(scene.masks.membrane, "all"));
verifyGreaterThan(testCase, scene.summary.max_abs_pressure, 0);

outDir = "C:\temp";
if ~isfolder(outDir)
    mkdir(outDir);
end
gifPath = string(tempname(outDir)) + ".gif";
testCase.addTeardown(@() deleteIfExists(gifPath));

info = writeDrumHighOrderImpedanceGif(scene, gifPath, "DelayTime", 0.01);

verifyEqual(testCase, info.kind, "drum_high_order_impedance_gif");
verifyEqual(testCase, info.num_frames, 4);
verifyTrue(testCase, isfile(gifPath));

frames = imfinfo(gifPath);
verifyEqual(testCase, numel(frames), 4);
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
