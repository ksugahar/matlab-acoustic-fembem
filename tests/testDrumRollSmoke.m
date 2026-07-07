function tests = testDrumRollSmoke
%TESTDRUMROLLSMOKE Fast drum-roll smoke for the tests/ lane.
%
% The heavy drum-roll validation -- directionality, full movies, and the
% drum+scatterer scene -- lives in validation_test/ (run before release).
% This keeps a quick check in the fast lane: the two-spot CQ drum roll still
% runs, alternates A/B, and is causal + real at minimal resolution.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDrumRollRunsCausalRealAtMinimalResolution(testCase)
result = drumRollConvolutionQuadrature( ...
    "NumBeats", 2, "TailTime", 1.0, "TimeStep", 0.5);   % default 2 listeners (directionality)

verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, cellstr(result.beatSpot(:)).', {'A', 'B'});   % two-spot alternation
verifyGreaterThan(testCase, result.summary.num_time, 3);
verifyGreaterThan(testCase, result.summary.max_abs_pressure, 0);

scale = result.summary.max_abs_pressure;
verifyLessThan(testCase, abs(result.pressure(1, 1)), 1e-3 * scale);                    % causal
verifyLessThan(testCase, result.summary.max_imag_pressure_before_real, 1e-8 * scale);  % real
verifyLessThan(testCase, result.summary.max_relative_residual, 1e-6);                  % well solved
end
