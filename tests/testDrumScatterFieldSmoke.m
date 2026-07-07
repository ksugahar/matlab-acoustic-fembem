function tests = testDrumScatterFieldSmoke
%TESTDRUMSCATTERFIELDSMOKE Fast MATLAB-only drum + sphere scatterer field smoke.
%
% The heavier drum + scatterer movie lives in validation_test/testDrumScatterField.
% This keeps a quick check in the fast lane: the two-body (cylinder drum + sphere)
% CQ field is built entirely in MATLAB (no Gmsh), the two bodies are auto-detected
% from the combined .vol, and the radiated + scattered pressure is finite / real /
% residual-small on a small grid.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDrumScatterFieldTwoBodyCausalReal(testCase)
field = drumScatterField("NumBeats", 2, "TimeStep", 0.4, "TailTime", 1.0, "NumGrid", 36);

verifyEqual(testCase, field.kind, "drum_scatter_two_body_field");
verifyEqual(testCase, field.status, "ok");
verifySize(testCase, field.pressure, [36, 36, field.summary.num_time]);

% two disjoint bodies detected; the sphere floats above and over the drum rim
verifyTrue(testCase, field.checks.two_bodies_detected);
verifyTrue(testCase, field.checks.scatterer_above_drum);
verifyLessThan(testCase, abs(field.scatterer_center(1) - field.drum_radius), 0.35*field.drum_radius);

% struck A / B alternating on the +x / -x drumhead spots
verifyEqual(testCase, cellstr(field.beatSpot(:)).', {'A', 'B'});
verifyTrue(testCase, field.checks.alternating_two_spot);

% inherited CQ health + physical field
verifyTrue(testCase, field.checks.finite_pressure);
verifyTrue(testCase, field.checks.real_time_response);
verifyTrue(testCase, field.checks.cq_residuals_small);
verifyGreaterThan(testCase, field.summary.max_abs_pressure, 0);

% both bodies masked out (NaN), exterior grid points finite
verifyTrue(testCase, any(field.mask_inside, "all"));
inflated = repmat(field.mask_inside, 1, 1, field.summary.num_time);
verifyTrue(testCase, all(isnan(field.pressure(inflated))));
end
