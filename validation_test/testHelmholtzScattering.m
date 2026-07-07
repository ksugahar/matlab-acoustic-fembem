function tests = testHelmholtzScattering
%TESTHELMHOLTZSCATTERING Acoustic Helmholtz single-layer rung, 3-way validated.
%
% Gates locked from the 2026-07-03 measurements: analytic references
% (interior point source = exact, sound-soft sphere partial-wave series),
% this repo's Galerkin Helmholtz path, and the committed ngsolve.bem
% references. All analytic deviations are faceted-geometry dominated
% (coarse -> fine improves ~x2.5, matching the O(h^2) faceting class); the
% two codes agree 10-30x tighter than either agrees with the true sphere.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "validation"));
end


function testPointSourceGateConvergesWithMesh(testCase)
% EXACT gate: boundary data from an interior point source must reproduce
% the same field at exterior points (uniqueness + radiation condition).
% measured: coarse k=2 max 7.65e-2, fine k=2 max 3.10e-2, fine k=0.5 1.63e-2
x0 = [0.3 0.2 -0.25];
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
allErr = struct();
for name = ["unit_sphere_coarse", "unit_sphere_fine"]
    surface = fixtureSurface(name + ".vol");
    for k = [0.5 2.0]
        g = acousticPointSource(k, x0, surface.vtx);
        sol = singleLayerDirichletSolve(surface, g, ...
            "Wavenumber", k, "QuadratureOrder", 7);
        verifyEqual(testCase, sol.status, "ok");
        verifyEqual(testCase, sol.kind, ...
            "helmholtz_single_layer_exterior_dirichlet_solve");
        verifyFalse(testCase, isreal(sol.q));
        pex = acousticPointSource(k, x0, probes);
        e = max(abs(sol.potentialAt(probes) - pex) ./ abs(pex));
        allErr.(name + "_k" + strrep(string(k), ".", "p")) = e;
    end
end
verifyLessThan(testCase, allErr.unit_sphere_fine_k0p5, 0.025);
verifyLessThan(testCase, allErr.unit_sphere_fine_k2, 0.05);
verifyLessThan(testCase, allErr.unit_sphere_fine_k0p5, ...
    allErr.unit_sphere_coarse_k0p5);
verifyLessThan(testCase, allErr.unit_sphere_fine_k2, ...
    allErr.unit_sphere_coarse_k2);
end


function testSoftSphereScatteringMatchesSeries(testCase)
% physics gate: plane wave on the sound-soft unit sphere vs the analytic
% partial-wave series. measured fine: k=0.5 1.37e-2, k=2.0 4.00e-2.
surface = fixtureSurface("unit_sphere_fine.vol");
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
bands = struct("k0p5", 0.025, "k2", 0.06);
for k = [0.5 2.0]
    g = -exp(1i * k * surface.vtx(:, 3));
    sol = singleLayerDirichletSolve(surface, g, ...
        "Wavenumber", k, "QuadratureOrder", 7);
    ref = softSphereScattering(k, 1.0, probes);
    verifyLessThan(testCase, ref.truncationTail, 1e-12);
    e = max(abs(sol.potentialAt(probes) - ref.scattered) ./ abs(ref.scattered));
    verifyLessThan(testCase, e, bands.("k" + strrep(string(k), ".", "p")));
end

% the series itself satisfies the soft boundary condition on the true sphere
theta = linspace(0.1, pi - 0.1, 7).';
onSphere = [sin(theta), zeros(size(theta)), cos(theta)];
ref = softSphereScattering(2.0, 1.0, onSphere);
verifyLessThan(testCase, max(abs(ref.total)), 1e-9);
end


function testOperatorSpectralLambda0(testCase)
% unit-sphere single-layer eigenvalue on the constant mode:
% lambda_0 = sin(k) e^{+ik} / k (the imaginary sign pins the e^{+ikr}
% convention). measured fine: k=0.5 rel 9.14e-3, k=2.0 rel 2.57e-2.
surface = fixtureSurface("unit_sphere_fine.vol");
n = size(surface.vtx, 1);
one = ones(n, 1);
[M, ~] = SurfaceP1Space(surface).mass();
bands = struct("k0p5", 0.02, "k2", 0.04);
for k = [0.5 2.0]
    op = GalerkinSingleLayer(surface, "Wavenumber", k, "QuadratureOrder", 7);
    lam0 = (one' * (op.matrix * one)) / (one' * (M * one));
    ana = sin(k) * exp(1i * k) / k;
    verifyLessThan(testCase, abs(lam0 - ana) / abs(ana), ...
        bands.("k" + strrep(string(k), ".", "p")));
    verifyGreaterThan(testCase, imag(lam0), 0);   % e^{+ikr}, not e^{-ikr}
end
end


function testLowFrequencyLimitMatchesLaplaceSolve(testCase)
% k -> 0 continuity of the solve AND the exterior evaluator (the expm1
% split makes the limit exact by construction).
% measured: density 7.1e-10, probe potential 5.3e-11.
surface = fixtureSurface("unit_sphere_coarse.vol");
x0 = [0.3 0.2 -0.25];
probes = [2 0 0; 0 0 3];
g = acousticPointSource(0, x0, surface.vtx);
s0 = singleLayerDirichletSolve(surface, g, "QuadratureOrder", 7);
se = singleLayerDirichletSolve(surface, g, ...
    "Wavenumber", 1e-9, "QuadratureOrder", 7);
verifyEqual(testCase, s0.kind, "laplace_single_layer_exterior_dirichlet_solve");
verifyTrue(testCase, s0.checks.densityTypeMatchesKernel);
verifyTrue(testCase, se.checks.densityTypeMatchesKernel);
verifyLessThan(testCase, norm(se.q - s0.q) / norm(s0.q), 1e-8);
verifyLessThan(testCase, ...
    max(abs(se.potentialAt(probes) - s0.potentialAt(probes))), 1e-9);
end


function testRigidScatteringAndIrregularFrequencies(testCase)
% sound-hard scattering via the total-field K equation
% (1/2 M - K_k) t = M g_inc, and the classic irregular-frequency lesson:
% at kR = pi (first interior Dirichlet eigenvalue) the equation is
% singular; the DISCRETE condition number stays benign (~29, the faceted
% eigenvalue shifts) while the solution is ~100% wrong - only the
% analytic gate sees it. CHIEF (interior null-field rows, least squares)
% restores regular-class accuracy. measured (fine, gss 7):
% k=2 direct 5.65e-2; k=pi direct 9.57e-1 -> chief 8.04e-2.
surface = fixtureSurface("unit_sphere_fine.vol");
probes = [2 0 0; 0 0 3; -1.2 1.6 0];

refReg = rigidSphereScattering(2.0, 1.0, probes);
verifyLessThan(testCase, refReg.truncationTail, 1e-12);
solReg = rigidScatteringSolve(surface, "Wavenumber", 2.0, ...
    "QuadratureOrder", 7);
verifyEqual(testCase, solReg.status, "ok");
verifyEqual(testCase, solReg.kind, ...
    "rigid_scattering_total_field_double_layer_solve");
eReg = max(abs(solReg.scatteredAt(probes) - refReg.scattered) ...
    ./ abs(refReg.scattered));
verifyLessThan(testCase, eReg, 0.08);

refIrr = rigidSphereScattering(pi, 1.0, probes);
solDir = rigidScatteringSolve(surface, "Wavenumber", pi, ...
    "QuadratureOrder", 7);
eDir = max(abs(solDir.scatteredAt(probes) - refIrr.scattered) ...
    ./ abs(refIrr.scattered));
verifyGreaterThan(testCase, eDir, 0.5);           % breakdown LOCKED
verifyLessThan(testCase, solDir.conditionNumber, 100);   % ...while cond benign
solChf = rigidScatteringSolve(surface, "Wavenumber", pi, ...
    "QuadratureOrder", 7, "Method", "chief");
eChf = max(abs(solChf.scatteredAt(probes) - refIrr.scattered) ...
    ./ abs(refIrr.scattered));
verifyLessThan(testCase, eChf, 0.12);             % CHIEF rescues
end


function testThreeWayNgsolveCrossCheck(testCase)
% standing regression from the committed ngbem Helmholtz .mat artifacts:
% operator, density, and probe-point agreement between the two codes must
% stay 10-30x tighter than the faceting-level analytic deviation.
% measured max over fixtures/cases: operator 8.7e-3, cross-code 6.2e-3.
for name = ["unit_sphere_coarse", "unit_sphere_fine"]
    report = verifyHelmholtzAgainstNgsolve(fixtureFile(name + ".vol"));
    verifyEqual(testCase, report.status, "ok");
    for c = 1:numel(report.cases)
        r = report.cases(c);
        verifyTrue(testCase, r.checks.conventionPinned);
        verifyTrue(testCase, r.checks.doubleLayerConventionPinned);
        verifyLessThan(testCase, r.operatorRelDiff, 2e-2);
        verifyLessThan(testCase, r.doubleLayerRelDiff, 1e-2);
        verifyLessThan(testCase, r.pointSourceProbeCrossCode, 2e-2);
        verifyLessThan(testCase, r.planeWaveProbeCrossCode, 2e-2);
        verifyLessThan(testCase, r.rigidTraceCrossCode, 1e-3);
        verifyLessThan(testCase, r.rigidProbeCrossCode, 1e-3);
        verifyLessThan(testCase, r.referenceIntorderConvergenceV, 1e-6);
    end
end
end


function surface = fixtureSurface(name)
mesh = VolMesh(fixtureFile(name));
surface = mesh.boundary();
end


function file = fixtureFile(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
file = fullfile(repoRoot, "fixtures", "mesh_topology", name);
end
