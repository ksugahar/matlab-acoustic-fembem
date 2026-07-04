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
