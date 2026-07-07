function tests = testVolTdBemConvolutionQuadrature
%TESTVOLTDBEMCONVOLUTIONQUADRATURE .vol boundary P1 CQ time-domain BEM.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testBdf1CqSingleLayerTdBemRunsOnVolBoundary(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol");

result = volTdBemConvolutionQuadrature(volFile, ...
    "NumTime", 8, ...
    "TimeStep", 0.6, ...
    "QuadratureOrder", 1, ...
    "Method", "BDF1");

verifyEqual(testCase, result.kind, "vol_p1_td_bem_convolution_quadrature_response");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.method, "lubich_bdf1_cq_laplace_domain_single_layer_bem");
verifyEqual(testCase, size(result.pressure, 1), 8);
verifyEqual(testCase, result.summary.num_cq_laplace_solves, 8);
verifyTrue(testCase, result.checks.vol_mesh_tri_tet);
verifyTrue(testCase, result.checks.p1_boundary_bem);
verifyTrue(testCase, result.checks.laplace_parameters_positive_real);
verifyTrue(testCase, result.checks.cq_residuals_small);
verifyTrue(testCase, result.checks.real_time_response);
verifyTrue(testCase, result.checks.not_frequency_sweep_ifft);
verifyGreaterThan(testCase, result.summary.max_abs_pressure, 0);
verifyLessThan(testCase, result.summary.max_relative_residual, 1e-10);
verifyLessThan(testCase, result.summary.max_imag_pressure_before_real, ...
    1e-8 * max(1, result.summary.max_abs_pressure));
end


function testBdf2CqUsesPositiveLaplaceParameters(testCase)
result = volTdBemConvolutionQuadrature("", ...
    "NumTime", 8, ...
    "TimeStep", 0.6, ...
    "QuadratureOrder", 1, ...
    "Method", "BDF2");

verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.method, "lubich_bdf2_cq_laplace_domain_single_layer_bem");
verifyTrue(testCase, all(real(result.cqLaplaceParameter) > 0));
verifyTrue(testCase, result.checks.cq_residuals_small);
verifyGreaterThan(testCase, result.summary.max_abs_density, 0);
end


function testBoundaryTimeDataSizeIsChecked(testCase)
verifyError(testCase, @() volTdBemConvolutionQuadrature("", ...
    "NumTime", 8, ...
    "BoundaryTimeData", zeros(8, 2)), ...
    "volTdBemConvolutionQuadrature:BoundaryTimeData");
end


function testCqSingleLayerRhsIsMassConsistent(testCase)
% Regression lock for the boundary-mass fix.  The CQ single layer solves the
% GALERKIN boundary integral equation  V(s) q = M ghat, so its OWN primitives
% (laplaceSingleLayerGalerkin + laplaceSingleLayerPotential, evaluated at the
% imaginary node s = -1i c0 k exactly as the solver does) must reproduce the
% analytic soft-sphere scattered field WITH the boundary P1 mass M, and be
% grossly wrong WITHOUT it.  Dropping M -- the pre-fix bug -- scaled the
% scattered amplitude ~12x the analytic value; the auto-scaled movie hid it,
% and every shape/causality/ratio check passed regardless.
%
% This is an OPERATOR-level analytic anchor by necessity: the exterior
% first-kind single layer is singular at the interior Dirichlet eigenvalues
% kR = n*pi, so a full time-domain response cannot be cleanly cross-checked by
% a real-frequency inverse FFT (any pulse with energy near kR = pi blows the
% reference up).  k = 1.8 < pi is below the first irregular frequency, giving a
% clean direct comparison to the partial-wave series on the movie's own mesh.
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_fine.vol");
surface = VolMesh(volFile).boundary();
[massB, ~] = SurfaceP1Space(surface).mass();

c0 = 1.0; k = 1.8; order = 7;
sNode = -1i * c0 * k;                             % laplaceSingleLayerGalerkin(-i c0 k) == Helmholtz(k)
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
V = laplaceSingleLayerGalerkin(surface, sNode, c0, order);
Spot = laplaceSingleLayerPotential(surface, probes, sNode, c0, order);
g = -exp(1i * k * surface.vtx(:, 3));            % soft-sphere Dirichlet trace
ref = softSphereScattering(k, 1.0, probes);
verifyLessThan(testCase, ref.truncationTail, 1e-12);

pWithMass    = Spot * (V \ (massB * g));
pWithoutMass = Spot * (V \ g);
relWithMass    = max(abs(pWithMass    - ref.scattered) ./ abs(ref.scattered));
relWithoutMass = max(abs(pWithoutMass - ref.scattered) ./ abs(ref.scattered));
verifyLessThan(testCase, relWithMass, 6e-2);     % measured 3.7e-2 (mass-consistent)
verifyGreaterThan(testCase, relWithoutMass, 1.0);% measured ~12 (the pre-fix bug)
end
