function tests = testVolFemBemCoupledConvolutionQuadrature
%TESTVOLFEMBEMCOUPLEDCONVOLUTIONQUADRATURE Volume P1 FEM + BEM CQ coupling.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testBdf1CoupledVolumeFemBemCqRunsOnInteriorNodeVol(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "four_tet_interior_node.vol");

result = volFemBemCoupledConvolutionQuadrature(volFile, ...
    "NumTime", 8, ...
    "TimeStep", 0.6, ...
    "QuadratureOrder", 1, ...
    "Method", "BDF1");

verifyEqual(testCase, result.kind, "vol_p1_volume_fem_boundary_bem_coupled_cq_response");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.method, "lubich_bdf1_cq_volume_p1_fem_johnson_nedelec_calderon_bem_coupling");
verifyEqual(testCase, result.couplingForm, "JohnsonNedelec");
verifyEqual(testCase, result.summary.num_coupled_cq_laplace_solves, 8);
verifyEqual(testCase, result.summary.num_interior_only_nodes, 1);
verifyTrue(testCase, result.checks.vol_mesh_tri_tet);
verifyTrue(testCase, result.checks.p1_volume_fem);
verifyTrue(testCase, result.checks.p1_boundary_bem);
verifyTrue(testCase, result.checks.johnson_nedelec_calderon_form);
verifyTrue(testCase, result.checks.double_layer_k_included);
verifyTrue(testCase, result.checks.laplace_parameters_positive_real);
verifyTrue(testCase, result.checks.coupled_residuals_small);
verifyTrue(testCase, result.checks.real_interior_response);
verifyTrue(testCase, result.checks.real_exterior_response);
verifyTrue(testCase, result.checks.not_exterior_only_td_bem);
verifyGreaterThan(testCase, result.summary.max_abs_interior_pressure, 0);
verifyGreaterThan(testCase, result.summary.max_abs_exterior_pressure, 0);
verifyGreaterThan(testCase, result.summary.max_double_layer_frobenius_norm, 0);
verifyLessThan(testCase, result.summary.max_relative_residual, 1e-10);
end


function testBdf2CoupledCqKeepsPositiveLaplaceParameters(testCase)
result = volFemBemCoupledConvolutionQuadrature("", ...
    "NumTime", 8, ...
    "TimeStep", 0.6, ...
    "QuadratureOrder", 1, ...
    "Method", "BDF2");

verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.method, "lubich_bdf2_cq_volume_p1_fem_johnson_nedelec_calderon_bem_coupling");
verifyTrue(testCase, all(real(result.cqLaplaceParameter) > 0));
verifyTrue(testCase, result.checks.double_layer_k_included);
verifyTrue(testCase, result.checks.coupled_residuals_small);
verifyGreaterThan(testCase, result.summary.max_abs_boundary_density, 0);
end


function testSingleLayerTeachingFormRemainsAvailableForRegression(testCase)
result = volFemBemCoupledConvolutionQuadrature("", ...
    "NumTime", 8, ...
    "TimeStep", 0.6, ...
    "QuadratureOrder", 1, ...
    "CouplingForm", "SingleLayerTeaching");

verifyEqual(testCase, result.couplingForm, "SingleLayerTeaching");
verifyEqual(testCase, result.method, "lubich_bdf1_cq_volume_p1_fem_single_layer_teaching_bem_coupling");
verifyFalse(testCase, result.checks.johnson_nedelec_calderon_form);
verifyFalse(testCase, result.checks.double_layer_k_included);
verifyEqual(testCase, result.summary.max_double_layer_frobenius_norm, 0);
verifyEqual(testCase, result.status, "needs_attention");
end


function testVolumeSourceTimeDataLengthIsChecked(testCase)
verifyError(testCase, @() volFemBemCoupledConvolutionQuadrature("", ...
    "NumTime", 8, ...
    "VolumeSourceTimeData", zeros(7, 1)), ...
    "volFemBemCoupledConvolutionQuadrature:VolumeSourceTimeData");
end
