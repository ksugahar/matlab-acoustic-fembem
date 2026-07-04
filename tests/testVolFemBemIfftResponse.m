function tests = testVolFemBemIfftResponse
%TESTVOLFEMBEMIFFTRESPONSE .vol P1 FEM/BEM frequency sweep to time response.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testUnitTetraFrequencyIfftResponseIsSolverReady(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol");

result = volFemBemIfftResponse(volFile, ...
    "NumTime", 16, ...
    "FinalTime", 24, ...
    "QuadratureOrder", 1);

verifyEqual(testCase, result.kind, "vol_p1_fem_bem_frequency_ifft_time_response");
verifyEqual(testCase, result.status, "ok");
verifyEqual(testCase, result.method, "frequency_domain_p1_fem_p1_bem_plus_inverse_fft");
verifyEqual(testCase, size(result.pressure, 1), 16);
verifyEqual(testCase, result.summary.num_frequency_solves, 8);
verifyTrue(testCase, result.checks.vol_mesh_tri_tet);
verifyTrue(testCase, result.checks.p1_volume_fem);
verifyTrue(testCase, result.checks.p1_boundary_bem);
verifyTrue(testCase, result.checks.frequency_solves_ok);
verifyTrue(testCase, result.checks.ifft_is_real);
verifyTrue(testCase, result.checks.nyquist_bin_zeroed);
verifyGreaterThan(testCase, result.summary.max_abs_pressure, 0);
verifyLessThan(testCase, result.summary.max_solve_residual, 1e-10);
end


function testOddTimeCountRejected(testCase)
verifyError(testCase, @() volFemBemIfftResponse("", "NumTime", 15), ...
    "volFemBemIfftResponse:NumTime");
end
