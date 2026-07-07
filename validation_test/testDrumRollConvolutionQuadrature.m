function tests = testDrumRollConvolutionQuadrature
%TESTDRUMROLLCONVOLUTIONQUADRATURE Two-spot alternating strikes (a drum roll) by CQ.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDrumRollAlternatesTwoSpotsCausalReal(testCase)
result = drumRollConvolutionQuadrature();

verifyEqual(testCase, result.kind, "drum_roll_two_spot_alternating_cq_time_response");
verifyEqual(testCase, result.status, "ok");

% three beats, struck A / B / A at the +x and -x drumhead spots (cylinder top face)
verifyEqual(testCase, result.summary.num_time, 20);
verifyEqual(testCase, cellstr(result.beatSpot(:)).', {'A', 'B', 'A'});
verifyTrue(testCase, result.checks.alternating_two_spot);
verifyGreaterThan(testCase, result.strikeSpotA(1), 0);   % +x on the drumhead
verifyLessThan(testCase, result.strikeSpotB(1), 0);      % -x on the drumhead
verifyEqual(testCase, result.strikeSpotA(3), result.strikeSpotB(3));  % same drumhead height
verifyTrue(testCase, result.checks.drumhead_strikes_on_top_face);     % a cylinder drum, not a sphere

% inherited CQ health: causal, real, well conditioned
verifyGreaterThan(testCase, result.summary.max_abs_pressure, 0);
verifyLessThan(testCase, result.summary.max_relative_residual, 1e-6);
scale = result.summary.max_abs_pressure;
verifyLessThan(testCase, abs(result.pressure(1, 1)), 1e-3 * scale);          % causal
verifyLessThan(testCase, result.summary.max_imag_pressure_before_real, 1e-8 * scale);  % real
end


function testDrumRollIsDirectionalOnTheAStruckSide(testCase)
result = drumRollConvolutionQuadrature();

% the odd-beat (A) taps are clearly loudest on the A (+x) listener side
d = result.directional;
verifyTrue(testCase, result.checks.directional_drum_roll);
verifyTrue(testCase, d.A_louder_on_A_side);
verifyGreaterThan(testCase, d.A_side_energy_ratio, 2.0);
verifyGreaterThan(testCase, d.A_energy_at_listenerA, d.A_energy_at_listenerB);
end


function testDrumRollHonoursCustomListeners(testCase)
obs = [3 0 0; 0 0 3];
result = drumRollConvolutionQuadrature( ...
    "NumBeats", 2, "ObservationPoints", obs);

verifyEqual(testCase, cellstr(result.beatSpot(:)).', {'A', 'B'});
verifyEqual(testCase, result.listeners, obs);
verifySize(testCase, result.pressure, [result.summary.num_time, 2]);
verifyEqual(testCase, result.status, "ok");
end
