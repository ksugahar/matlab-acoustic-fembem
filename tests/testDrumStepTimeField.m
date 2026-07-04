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
verifyEqual(testCase, scene.geometry.radiation_model, "one-sided baffled top-head radiation");
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
verifyTrue(testCase, info.flip_vertical);
verifyTrue(testCase, isfile(gifPath));

frames = imfinfo(gifPath);
verifyEqual(testCase, numel(frames), 4);
end


function testReducedFemBemDrumHasBottomRadiation(testCase)
scene = drumFemBemCoupledDemo( ...
    "NumX", 42, ...
    "NumZ", 42, ...
    "NumTime", 14, ...
    "TMax", 1.0e-3, ...
    "NumSourceRadial", 4, ...
    "NumSourceAzimuth", 8, ...
    "NumSideAxial", 4);

verifyEqual(testCase, scene.kind, "drum_reduced_fem_bem_coupled_scene");
verifyEqual(testCase, scene.status, "ok");
verifyTrue(testCase, contains(scene.coupling.kind, "fem"));
verifyTrue(testCase, contains(scene.boundary_type, "high-order impedance"));
verifyTrue(testCase, contains(scene.bem.observation_rule, "every exterior air observation point"));
verifyTrue(testCase, scene.checks.lower_half_wave_present);
verifyTrue(testCase, scene.checks.internal_cavity_coupled);
verifyTrue(testCase, scene.checks.top_source_reaches_lateral_exterior);
verifyTrue(testCase, scene.checks.top_source_reaches_lower_exterior);
verifyTrue(testCase, scene.checks.bottom_source_reaches_upper_exterior);
verifyTrue(testCase, scene.checks.side_source_reaches_upper_exterior);
verifyTrue(testCase, scene.checks.side_source_reaches_lower_exterior);
verifyTrue(testCase, any(scene.masks.bottom_membrane, "all"));
verifyTrue(testCase, any(scene.masks.cavity, "all"));
verifyGreaterThan(testCase, scene.summary.bottom_peak_acceleration, 0);
verifyGreaterThan(testCase, scene.summary.cavity_peak_pressure, 0);
verifyGreaterThan(testCase, scene.summary.bem_cross_direction.top_to_lateral_exterior, 0);
verifyGreaterThan(testCase, scene.summary.bem_cross_direction.top_to_lower_exterior, 0);
verifyGreaterThan(testCase, scene.summary.bem_cross_direction.side_to_upper_exterior, 0);

outDir = "C:\temp";
if ~isfolder(outDir)
    mkdir(outDir);
end
gifPath = string(tempname(outDir)) + ".gif";
testCase.addTeardown(@() deleteIfExists(gifPath));

info = writeDrumHighOrderImpedanceGif(scene, gifPath, ...
    "DelayTime", 0.01, ...
    "OutputSize", [64, 64]);

verifyEqual(testCase, info.kind, "drum_high_order_impedance_gif");
verifyTrue(testCase, info.axis_equal);
verifyTrue(testCase, info.flip_vertical);
verifyEqual(testCase, info.output_size, [64, 64]);
verifyTrue(testCase, isfile(gifPath));

frames = imfinfo(gifPath);
verifyEqual(testCase, double(frames(1).Height), 64);
verifyEqual(testCase, double(frames(1).Width), 64);
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
