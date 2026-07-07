function tests = testFemBemHelmholtzCoupling
%TESTFEMBEMHELMHOLTZCOUPLING Pure-MATLAB acoustic FEM/BEM coupled solve.
%
% The three coupled acoustic cases (2026-07-03 measurements), all solved
% by femBemCoupledSolve alone - interior P1 Helmholtz FEM + exterior P1
% Galerkin BEM, no ABC anywhere (the BEM row IS the radiation condition):
%
%   case 1  k -> 0 regression against the verified Laplace coupled solve
%   case 2  acoustic invisibility (k1 = k0, rho1 = rho0): the EXACT null
%           gate - interior == incident plane wave, scattered == 0, up to
%           discretization, converging under mesh refinement
%   case 3  Anderson (1950) fluid-sphere transmission (c1/c0 = 0.7,
%           rho1/rho0 = 1.2) against the partial-wave series, converging
%           under mesh refinement
%
% Plus the sphere spectral gates of the new Helmholtz double layer K_k.
% Fixtures: unit_sphere_fine (maxh 0.3) and unit_ball_maxh018 - the pair
% turns the P1 (k h)^2 resolution limit into a measured convergence
% assertion instead of a hidden loose band.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testHelmholtzDoubleLayerSphereSpectralGates(testCase)
% K_k[Y_l] = 1/2 + 1i k^2 j_l(k) h_l'(k) on the unit sphere (our
% outward-normal PV convention; Laplace limit -1/(2(2l+1))).
% measured gss7 fine: k=0.5 -> 3.3e-3 / 4.7e-3; k=2 -> 2.9e-2 / 2.3e-2;
% k->0 vs Laplace K: 5.4e-28.
surface = fixtureSurface("unit_sphere_fine.vol");
n = size(surface.vtx, 1);
one = ones(n, 1);
zc = surface.vtx(:, 3);
[M, ~] = SurfaceP1Space(surface).mass();

sphJ = @(l, x) sqrt(pi ./ (2 * x)) .* besselj(l + 0.5, x);
sphH = @(l, x) sqrt(pi ./ (2 * x)) .* (besselj(l + 0.5, x) + 1i * bessely(l + 0.5, x));
dH = @(l, x) sphH(l - 1, x) - (l + 1) ./ x .* sphH(l, x);

bands = struct("k0p5", 0.01, "k2", 0.05);
for k = [0.5 2.0]
    K = GalerkinDoubleLayer(surface, "Wavenumber", k, "QuadratureOrder", 7);
    band = bands.("k" + strrep(string(k), ".", "p"));
    lam0 = (one' * (K.matrix * one)) / (one' * (M * one));
    ana0 = 0.5 + 1i * k^2 * sphJ(0, k) * dH(0, k);
    verifyLessThan(testCase, abs(lam0 - ana0) / abs(ana0), band);
    lam1 = (zc' * (K.matrix * zc)) / (zc' * (M * zc));
    ana1 = 0.5 + 1i * k^2 * sphJ(1, k) * dH(1, k);
    verifyLessThan(testCase, abs(lam1 - ana1) / abs(ana1), band);
end

K0 = GalerkinDoubleLayer(surface, "QuadratureOrder", 7);
Ke = GalerkinDoubleLayer(surface, "Wavenumber", 1e-9, "QuadratureOrder", 7);
verifyLessThan(testCase, ...
    norm(Ke.matrix - K0.matrix, "fro") / norm(K0.matrix, "fro"), 1e-12);
end


function testCase1LowFrequencyLimitMatchesLaplace(testCase)
% measured: u 9.5e-10, lambda 1.8e-12, exterior 3.2e-10.
m = fixtureModel("unit_sphere_fine.vol");
probes = [2 0 0; 0 0 3];
s0 = femBemCoupledSolve(m);
sE = femBemCoupledSolve(m, "Wavenumber", 1e-9);
verifyEqual(testCase, s0.kind, "johnson_nedelec_coupled_fem_bem_solve");
verifyEqual(testCase, sE.kind, "johnson_nedelec_coupled_fem_bem_helmholtz_solve");
verifyEqual(testCase, s0.status, "ok");
verifyEqual(testCase, sE.status, "ok");
verifyLessThan(testCase, norm(sE.u - s0.u) / norm(s0.u), 1e-8);
verifyLessThan(testCase, norm(sE.lambda - s0.lambda) / norm(s0.lambda), 1e-10);
verifyLessThan(testCase, ...
    max(abs(sE.exteriorPotentialAt(probes) - s0.exteriorPotentialAt(probes))), 1e-8);
end


function testCase2AcousticInvisibilityExactNull(testCase)
% k1 = k0, rho1 = rho0: the sphere is acoustically invisible. The series
% itself is exact (A_l = 0 analytically after the log-derivative
% elimination; self-check 4.7e-8 incl. far probes); the coupled solve
% reproduces the null up to discretization and CONVERGES:
% fine 4.15e-2 / 2.57e-2  ->  finer ball 1.26e-2 / 7.83e-3.
k = 2.0;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];

mF = fixtureModel("unit_sphere_fine.vol");
ref = fluidSphereScattering(k, 1.0, [mF.mesh.vtx; probes]);
verifyLessThan(testCase, max(abs(ref.total - ref.incident)), 1e-6);

err = struct();
for name = ["unit_sphere_fine", "unit_ball_maxh018"]
    m = fixtureModel(name + ".vol");
    sol = femBemCoupledSolve(m, "Wavenumber", k, ...
        "VolumeSource", 0, "IncidentAmplitude", 1);
    verifyEqual(testCase, sol.status, "ok");
    verifyTrue(testCase, sol.checks.solutionFieldTypeMatchesKernel);
    pinc = exp(1i * k * m.mesh.vtx(:, 3));
    err.(name) = struct( ...
        "interior", norm(sol.u - pinc) / norm(pinc), ...
        "scattered", max(abs(sol.exteriorPotentialAt(probes))));
end
verifyLessThan(testCase, err.unit_sphere_fine.interior, 0.07);
verifyLessThan(testCase, err.unit_sphere_fine.scattered, 0.05);
verifyLessThan(testCase, err.unit_ball_maxh018.interior, 0.025);
verifyLessThan(testCase, err.unit_ball_maxh018.scattered, 0.015);
verifyLessThan(testCase, err.unit_ball_maxh018.interior, ...
    err.unit_sphere_fine.interior);
verifyLessThan(testCase, err.unit_ball_maxh018.scattered, ...
    err.unit_sphere_fine.scattered);
end


function testCase3AndersonFluidSphereTransmission(testCase)
% c1/c0 = 0.7 (k1 = k0/0.7), rho1/rho0 = 1.2 at k0 = 2 vs the Anderson
% partial-wave series. measured: interior/trace/probes(max)
% fine 1.30e-1 / 1.28e-1 / 2.18e-1  ->  finer 4.36e-2 / 4.24e-2 / 7.33e-2
% (the P1 (k1 h)^2 resolution class, converging ~h^2).
k0 = 2.0;
k1 = k0 / 0.7;
rhor = 1.2;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
refP = fluidSphereScattering(k0, 1.0, probes, ...
    "InteriorWavenumber", k1, "DensityRatio", rhor);
verifyLessThan(testCase, refP.truncationTail, 1e-6);
scatRef = refP.total - refP.incident;

err = struct();
for name = ["unit_sphere_fine", "unit_ball_maxh018"]
    m = fixtureModel(name + ".vol");
    sol = femBemCoupledSolve(m, "Wavenumber", k0, "InteriorWavenumber", k1, ...
        "DensityRatio", rhor, "VolumeSource", 0, "IncidentAmplitude", 1);
    verifyEqual(testCase, sol.status, "ok");
    ref = fluidSphereScattering(k0, 1.0, m.mesh.vtx, ...
        "InteriorWavenumber", k1, "DensityRatio", rhor);
    us = sol.exteriorPotentialAt(probes);
    err.(name) = struct( ...
        "interior", norm(sol.u - ref.total) / norm(ref.total), ...
        "probes", max(abs(us - scatRef) ./ abs(scatRef)));
end
verifyLessThan(testCase, err.unit_sphere_fine.interior, 0.20);
verifyLessThan(testCase, err.unit_ball_maxh018.interior, 0.07);
verifyLessThan(testCase, err.unit_ball_maxh018.probes, 0.12);
verifyLessThan(testCase, err.unit_ball_maxh018.interior, ...
    err.unit_sphere_fine.interior);
verifyLessThan(testCase, err.unit_ball_maxh018.probes, ...
    err.unit_sphere_fine.probes);
end


function surface = fixtureSurface(name)
surface = fixtureModel(name).surface;
end


function m = fixtureModel(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
m = FemBemModel(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end
