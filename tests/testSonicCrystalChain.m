function tests = testSonicCrystalChain
%TESTSONICCRYSTALCHAIN Multi-body scattering rung on the sphere chain.
%
% Five sound-soft spheres (R = 0.3, lattice constant d = 1.5) in free
% space - the sonic-crystal teaching geometry this first-order lane
% supports. Locks BOTH the multi-body machinery (4-leg cross-check via
% verifySonicCrystalChain) and the honest PHYSICS FINDING: a sparse
% free-space chain shows broadband sub-wavelength attenuation with NO
% Bragg stop band at k d = pi (both the BEM and the Foldy model agree) -
% the band gap of the COMSOL "Sonic Crystal" class model needs duct
% confinement / transverse periodicity, which is the declared next rung
% (Bloch unit cell + duct transmission FEM).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "validation"));
end


function testMultiBodyPointSourceExactGate(testCase)
% interior point source inside the MIDDLE sphere must be reproduced at
% exterior points by the 5-body solve (uniqueness gate; gss 3).
% measured 2026-07-03: k=2.0944 max 1.98e-2.
surface = chainSurface();
x0 = [0 0 0.05];
probes = [0 0 4.2; 0 0 5.0; 0.6 0 4.5];
for k = [1.0 2.0944]
    g = acousticPointSource(k, x0, surface.vtx);
    sol = singleLayerDirichletSolve(surface, g, ...
        "Wavenumber", k, "QuadratureOrder", 3);
    verifyEqual(testCase, sol.status, "ok");
    pex = acousticPointSource(k, x0, probes);
    e = max(abs(sol.potentialAt(probes) - pex) ./ abs(pex));
    verifyLessThan(testCase, e, 0.04);
end
end


function testFoldyAgreesAtLowFrequencyAndDegradesHonestly(testCase)
% the monopole-only Foldy reference is a LOW-k reference: chain total
% field agrees to a few percent in the sub-wavelength regime and the
% deviation grows monotonically with k as the neglected l >= 1
% single-sphere terms turn on. measured: 1.4e-2 (k=0.6) -> 2.3e-1 (k=3).
surface = chainSurface();
probes = [0 0 4.2; 0 0 5.0; 0.6 0 4.5];
centers = [zeros(5, 2), (-3:1.5:3).'];
devs = zeros(1, 3);
ks = [0.6 1.0 3.0];
for m = 1:3
    k = ks(m);
    sol = singleLayerDirichletSolve(surface, -exp(1i * k * surface.vtx(:, 3)), ...
        "Wavenumber", k, "QuadratureOrder", 3);
    ptot = sol.potentialAt(probes) + exp(1i * k * probes(:, 3));
    foldy = foldyPointScattering(k, 0.3, centers, probes);
    devs(m) = max(abs(ptot - foldy.total) ./ abs(foldy.total));
end
verifyLessThan(testCase, devs(1), 0.05);
verifyLessThan(testCase, devs(2), 0.06);
verifyGreaterThan(testCase, devs(3), devs(1));   % honest degradation lock
end


function testNoBraggStopBandInFreeSpaceChain(testCase)
% NEGATIVE-RESULT LOCK: the free-space chain does NOT develop a stop band
% at the Bragg wavenumber pi/d = 2.0944. Insertion loss stays a broadband
% sub-wavelength-attenuation plateau (measured 3.47..3.94 dB over
% k = 0.6..3.0, spread 0.47 dB). If this test ever fails toward a DEEP
% dip, the geometry/solver changed - not the physics of sparse free-space
% chains.
surface = chainSurface();
probes = [0 0 4.2; 0 0 5.0; 0.6 0 4.5];
ks = [0.6 1.4 2.2 3.0];
il = zeros(size(ks));
for m = 1:numel(ks)
    k = ks(m);
    sol = singleLayerDirichletSolve(surface, -exp(1i * k * surface.vtx(:, 3)), ...
        "Wavenumber", k, "QuadratureOrder", 3);
    ptot = sol.potentialAt(probes) + exp(1i * k * probes(:, 3));
    il(m) = -20 * log10(mean(abs(ptot)));
end
verifyGreaterThan(testCase, min(il), 2.5);   % attenuation IS present
verifyLessThan(testCase, max(il), 5.0);
verifyLessThan(testCase, max(il) - min(il), 1.5);   % ... and featureless
end


function testFourLegCrossCheck(testCase)
% standing regression from the committed NGSolve chain .mat: exact gate,
% ngbem cross-code, FEM methods-diverse leg, Foldy behavior, convention.
report = verifySonicCrystalChain();
verifyEqual(testCase, report.status, "ok");
for c = 1:numel(report.cases)
    r = report.cases(c);
    verifyTrue(testCase, r.checks.conventionPinned);
    verifyLessThan(testCase, r.pointSourceExactGate, 4e-2);
    verifyLessThan(testCase, r.planeWaveProbeCrossCode, 3e-2);
    verifyLessThan(testCase, r.planeWaveProbeVsFem, 0.10);
    verifyLessThan(testCase, r.referenceIntorderConvergenceV, 1e-6);
end
end


function surface = chainSurface()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
mesh = VolMesh(fullfile(repoRoot, "fixtures", "mesh_topology", ...
    "soft_sphere_chain_5.vol"));
surface = mesh.boundary();
end
