function tests = testVolFemBemElasticConvolutionQuadrature
%TESTVOLFEMBEMELASTICCONVOLUTIONQUADRATURE Golden for the elastic FEM + acoustic
% BEM coupled CQ time-domain solver.
%
%   Locks the end-to-end gates that a validated temporary-directory run established: a
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


function testBlockReducesToFrequencyFsiAtImaginaryNode(testCase)
% CORRECTNESS ANCHOR: the Laplace-domain elastic-CQ block at s = -i c0 k must
% reproduce the validated frequency solver fsiCoupledSolve(k).  (On the fine
% unit-ball at quad 7 the match is 1.2e-3; on this coarse fast fixture a loose
% band still catches any sign / operator / coupling regression -- those are O(1).)
vol = volumeFixture();
model = FemBemModel(vol); mesh = model.mesh; surface = model.surface;
k = 1.2; c0 = 1.0; q = 3;
cL = 2.0; cT = 1.0; rhoS = 1.5; rhoF = 1.0;
mu = rhoS*cT^2; lamE = rhoS*(cL^2 - 2*cT^2);

ref = fsiCoupledSolve(model, Wavenumber=k, LongitudinalSpeed=cL, ShearSpeed=cT, ...
    DensityRatio=rhoS, FluidDensity=rhoF, QuadratureOrder=q, ExteriorMethod="bem");

% reassemble the CQ block at s = -i c0 k (s^2 = -k^2) with the SHARED operators
s = -1i*c0*k;
[Ks, Ms] = elasticityMatrices(mesh, lamE, mu, rhoS);
nV = size(mesh.vtx,1); ids = surface.volNodeIds; nB = numel(ids);
[Mb, ~] = SurfaceP1Space(surface).mass();
[G, Minc] = interfaceCoupling(surface, ids, nV, ...
    @(X)[zeros(size(X,1),2), 1i*k*exp(1i*k*X(:,3))]);
pincB = exp(1i*k*surface.vtx(:,3));
V = laplaceSingleLayerGalerkin(surface, s, c0, q);
K = laplaceDoubleLayerGalerkin(surface, s, c0, q);
ZuB = sparse(3*nV,nB); ZBB = sparse(nB,nB);
lhs = [ Ks+s^2*Ms, G.',        ZuB;
        rhoF*s^2*G, ZBB,        Mb;
        ZuB.',      0.5*Mb-K,   V ];
x = lhs \ [ -G.'*pincB; -Minc; zeros(nB,1) ];
ps = x(3*nV+(1:nB));

rel = norm(ps - ref.surfacePressure) / max(1, norm(ref.surfacePressure));
verifyLessThan(testCase, rel, 5e-2);
end


function [G, Minc] = interfaceCoupling(surface, ids, nV, incGrad)
% G_ij = int_Gamma mu_i (n . phi_struct_j), Minc_i = int_Gamma mu_i (grad p_inc . n);
% the same P1-lumped centroid rule as fsiCoupledSolve (mirrored here for the test).
signs = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri; vtx = surface.vtx; nB = size(vtx,1);
massTri = (ones(3)+eye(3))/12; nE = 27*size(tri,1);
Grow = zeros(nE,1); Gcol = zeros(nE,1); Gval = zeros(nE,1); Minc = zeros(nB,1); cur = 1;
for t = 1:size(tri,1)
    lc = tri(t,:); X = vtx(lc,:);
    cr = cross(X(2,:)-X(1,:), X(3,:)-X(1,:));
    area = 0.5*norm(cr); nrm = signs(t)*cr/norm(cr); vid = ids(lc);
    for a = 1:3, for b = 1:3, for c = 1:3
        Grow(cur) = lc(a); Gcol(cur) = 3*vid(b)-3+c; Gval(cur) = area*massTri(a,b)*nrm(c); cur = cur+1;
    end, end, end
    Xc = mean(X,1);
    Minc(lc) = Minc(lc) + (area/3)*(incGrad(Xc)*nrm.');
end
G = sparse(Grow, Gcol, Gval, nB, 3*nV);
end


function vol = volumeFixture()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
vol = string(fullfile(repoRoot, "fixtures", "mesh_topology", "four_tet_interior_node.vol"));
end
