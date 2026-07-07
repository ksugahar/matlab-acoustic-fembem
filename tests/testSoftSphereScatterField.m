function tests = testSoftSphereScatterField
%TESTSOFTSPHERESCATTERFIELD Time-domain soft-sphere pulse scattering movie + anchor.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testScatteredAmplitudeMatchesAnalyticSeries(testCase)
% Amplitude anchor for the movie: the CQ solver solves, at every Laplace node,
% a frequency soft-sphere BEM with this single-layer operator (anchor test:
% CQ reduces to it at s = -i c0 k).  Here we confirm that operator reproduces
% the analytic partial-wave series in the movie's mid-band (k ~ 1.8), on the
% SAME mesh + probes the movie uses.  (testHelmholtzScattering covers k=0.5, 2.0;
% this fills the pulse-peak band the scatter movie actually contains.)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_fine.vol"));
surface = VolMesh(volFile).boundary();
probes = [2 0 0; 0 0 3; -1.2 1.6 0];

k = 1.8;
g = -exp(1i * k * surface.vtx(:, 3));
sol = singleLayerDirichletSolve(surface, g, "Wavenumber", k, "QuadratureOrder", 7);
ref = softSphereScattering(k, 1.0, probes);
verifyLessThan(testCase, ref.truncationTail, 1e-12);
relerr = max(abs(sol.potentialAt(probes) - ref.scattered) ./ abs(ref.scattered));
verifyLessThan(testCase, relerr, 6e-2);   % measured 3.6e-2
end


function testScatterFieldCausalRealBoundedArrival(testCase)
field = softSphereScatterField("NumTime", 22, "TimeStep", 0.32, "NumGrid", 40);

verifyEqual(testCase, field.kind, "soft_sphere_scatter_field");
verifyEqual(testCase, field.status, "ok");
verifySize(testCase, field.pressure, [40, 40, 22]);

% robust physics gates
verifyTrue(testCase, field.checks.finite_pressure);
verifyTrue(testCase, field.checks.real_time_response);
verifyTrue(testCase, field.checks.cq_residuals_small);

% the scattered back-reflection peaks within a few steps of the ray arrival,
% never before it (causality)
verifyTrue(testCase, field.checks.arrival_peak_causal);
verifyTrue(testCase, field.checks.arrival_peak_near_geometry);
verifyGreaterThanOrEqual(testCase, field.summary.measured_peak, ...
    field.summary.geometric_arrival - field.time_step);
verifyLessThanOrEqual(testCase, field.summary.measured_peak, ...
    field.summary.geometric_arrival + 3 * field.time_step);
verifyGreaterThan(testCase, field.summary.max_abs_pressure, 0);

% inside-the-sphere grid points are masked out (NaN), exterior points are finite
verifyTrue(testCase, any(field.mask_inside, "all"));
verifyTrue(testCase, all(isnan(field.pressure(repmat(field.mask_inside, 1, 1, 22)))));
end


function testScatterGifWritesAnimation(testCase)
field = softSphereScatterField("NumTime", 16, "TimeStep", 0.4, "NumGrid", 36);

outDir = tempdir;
if ~isfolder(outDir)
    mkdir(outDir);
end
gifPath = string(tempname(outDir)) + ".gif";
testCase.addTeardown(@() deleteIfExists(gifPath));

info = writeSoftSphereScatterGif(field, gifPath, ...
    "TimeIndices", [2, 6, 10, 14], "DelayTime", 0.01);

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
