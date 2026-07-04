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
