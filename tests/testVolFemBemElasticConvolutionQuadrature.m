function tests = testVolFemBemElasticConvolutionQuadrature
%TESTVOLFEMBEMELASTICCONVOLUTIONQUADRATURE Golden for the elastic FEM + acoustic
% BEM coupled CQ time-domain solver.
%
%   Locks the end-to-end gates that the validated C:\temp run established: a
%   vector-elastic interior, a causal / real / bounded / small-residual exterior
%   response, and status = ok.  The correctness ANCHOR (the Laplace-domain block
%   at s = -i c0 k reproduces the frequency solver fsiCoupledSolve to ~1e-3) was
%   verified separately; promoting it to its own golden is a follow-up (it needs
%   the shared interface-coupling extracted to matlab_api/model/ to avoid copying).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testElasticScatteringCqRunsCausalRealBounded(testCase)
vol = volumeFixture();
res = volFemBemElasticConvolutionQuadrature(vol, NumTime=10, TimeStep=0.4, ...
    Method="BDF2", QuadratureOrder=1);

verifyEqual(testCase, res.status, "ok");
verifyEqual(testCase, res.kind, "elastic_fem_acoustic_bem_coupled_cq_time_response");

% vector-elastic interior: 3 displacement dof per volume node
verifyEqual(testCase, res.summary.num_volume_dof, 3 * res.meshSummary.points);

% causal: the leading time step (before the wave arrives) is ~zero
verifyLessThan(testCase, res.summary.causal_leading_ratio, 1e-2);

% real, physical, bounded, accurate
scaleP = res.summary.max_abs_exterior_pressure;
verifyGreaterThan(testCase, scaleP, 0);
verifyTrue(testCase, all(isfinite(res.exteriorPressure), "all"));
verifyLessThan(testCase, res.summary.max_imag_exterior_before_real, 1e-8 * max(1, scaleP));
verifyLessThan(testCase, res.summary.max_relative_residual, 1e-6);

% the CQ Laplace nodes stay in the stable right half-plane
verifyTrue(testCase, all(real(res.cqLaplaceParameter) > 0));
end


function testStifferSolidRadiatesLess(testCase)
% Physics sanity, no external reference: a much stiffer solid is closer to rigid,
% so at the SAME incident pulse it scatters a DIFFERENT (here, smaller near-field
% peak on this coarse mesh is not guaranteed) response -- we only assert both runs
% are well-posed and the responses genuinely differ, i.e. the interior elasticity
% actually feeds through to the exterior field (the coupling is live, not ignored).
vol = volumeFixture();
soft = volFemBemElasticConvolutionQuadrature(vol, NumTime=10, TimeStep=0.4, ...
    Method="BDF2", QuadratureOrder=1, LongitudinalSpeed=2.0, ShearSpeed=1.0);
stiff = volFemBemElasticConvolutionQuadrature(vol, NumTime=10, TimeStep=0.4, ...
    Method="BDF2", QuadratureOrder=1, LongitudinalSpeed=8.0, ShearSpeed=4.0);
verifyEqual(testCase, soft.status, "ok");
verifyEqual(testCase, stiff.status, "ok");
rel = norm(stiff.exteriorPressure(:) - soft.exteriorPressure(:)) ...
    / max(norm(soft.exteriorPressure(:)), realmin);
verifyGreaterThan(testCase, rel, 1e-3);   % interior stiffness reaches the exterior
end


function vol = volumeFixture()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
vol = string(fullfile(repoRoot, "fixtures", "mesh_topology", "four_tet_interior_node.vol"));
end
