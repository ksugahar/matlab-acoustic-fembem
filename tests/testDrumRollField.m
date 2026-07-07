function tests = testDrumRollField
%TESTDRUMROLLFIELD Spatial x-z movie field of the two-spot drum roll.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testDrumRollFieldRadiatesTwoSpotAlternatingCausalReal(testCase)
field = drumRollField("NumBeats", 3, "TimeStep", 0.3, "NumGrid", 40);

verifyEqual(testCase, field.kind, "drum_roll_two_spot_field");
verifyEqual(testCase, field.status, "ok");
verifySize(testCase, field.pressure, [40, 40, field.summary.num_time]);

% struck A / B / A at the +x and -x poles
verifyEqual(testCase, cellstr(field.beatSpot(:)).', {'A', 'B', 'A'});
verifyTrue(testCase, field.checks.alternating_two_spot);
verifyGreaterThan(testCase, field.strikeSpotA(1), 0);        % +x pole
verifyLessThan(testCase, field.strikeSpotB(1), 0);           % -x pole

% inherited CQ health + physical field
verifyTrue(testCase, field.checks.finite_pressure);
verifyTrue(testCase, field.checks.real_time_response);
verifyTrue(testCase, field.checks.cq_residuals_small);
verifyGreaterThan(testCase, field.summary.max_abs_pressure, 0);

% inside the sphere is masked out (NaN), exterior grid points are finite
verifyTrue(testCase, any(field.mask_inside, "all"));
inflated = repmat(field.mask_inside, 1, 1, field.summary.num_time);
verifyTrue(testCase, all(isnan(field.pressure(inflated))));
end


function testDrumRollFieldGifWritesAnimation(testCase)
field = drumRollField("NumBeats", 2, "TimeStep", 0.4, "NumGrid", 36);

outDir = tempdir;
if ~isfolder(outDir)
    mkdir(outDir);
end
gifPath = string(tempname(outDir)) + ".gif";
testCase.addTeardown(@() deleteIfExists(gifPath));

info = writeSoftSphereScatterGif(field, gifPath, ...
    "TimeIndices", [2, 5, 8, 11], "DelayTime", 0.01);

verifyEqual(testCase, info.kind, "soft_sphere_scatter_gif");
verifyEqual(testCase, info.num_frames, 4);
verifyTrue(testCase, isfile(gifPath));

frames = imfinfo(gifPath);
verifyEqual(testCase, numel(frames), 4);
verifyEqual(testCase, double(frames(1).Width), 36);
verifyEqual(testCase, double(frames(1).Height), 36);
end


function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
