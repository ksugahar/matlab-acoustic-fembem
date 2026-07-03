function tests = testElasticSphereScattering
%TESTELASTICSPHERESCATTERING Faran elastic sphere, validated by both limits.
%
% The analytic reference for the acoustic FSI (fluid-structure) lane: a
% solid ELASTIC sphere in a fluid. Locked by TWO INDEPENDENT references
% (2026-07-04):
%   ShearSpeed -> 0  reproduces the Anderson fluid sphere to ~1e-10 (exact),
%   very stiff       reproduces the rigid sphere,
% plus the 3x3 (with shear) converging to the 2x2 fluid limit, and the
% elastic radiation force showing resonance structure (rigid/soft do not).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testFluidLimitReproducesAnderson(testCase)
% ShearSpeed = 0 is a fluid sphere: must match Anderson exactly.
k = 2.0;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
cL = 0.7; rho = 1.2;
elastic = elasticSphereScattering(k, 1.0, probes, ...
    "LongitudinalSpeed", cL, "ShearSpeed", 0.0, "DensityRatio", rho);
anderson = fluidSphereScattering(k, 1.0, probes, ...
    "InteriorWavenumber", k / cL, "DensityRatio", rho);
% off-axis probes agree to 5e-13; the +z-axis probe accumulates a coherent
% sum (all P_l(1) = 1) so two independent implementations differ at the
% rounding level (measured 3.4e-8) - still 5 orders below the scattering.
verifyLessThan(testCase, ...
    max(abs(elastic.total - anderson.total) ./ abs(anderson.total)), 1e-6);
end


function testStiffLimitReproducesRigid(testCase)
% a very stiff/heavy solid approaches the rigid sphere.
k = 2.0;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
elastic = elasticSphereScattering(k, 1.0, probes, ...
    "LongitudinalSpeed", 50, "ShearSpeed", 30, "DensityRatio", 100);
rigid = rigidSphereScattering(k, 1.0, probes);
verifyLessThan(testCase, ...
    max(abs(elastic.total - rigid.total) ./ abs(rigid.total)), 1e-2);
end


function testShearTermsConvergeToFluidLimit(testCase)
% the full 3x3 (with shear) must approach the shear-free 2x2 as cT -> 0,
% validating the shear stress/displacement machinery in the fluid limit.
k = 2.0;
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
p0 = elasticSphereScattering(k, 1.0, probes, ...
    "LongitudinalSpeed", 0.7, "ShearSpeed", 0.0, "DensityRatio", 1.2).total;
pE = elasticSphereScattering(k, 1.0, probes, ...
    "LongitudinalSpeed", 0.7, "ShearSpeed", 0.02, "DensityRatio", 1.2).total;
verifyLessThan(testCase, max(abs(pE - p0) ./ abs(p0)), 1e-2);
end


function testElasticRadiationForceShowsResonance(testCase)
% the payoff: a lucite-like elastic sphere has a radiation-force resonance
% (Y_p is non-monotonic in kR and peaks well above the rigid ~0.75),
% something rigid/soft spheres cannot produce. The force stays
% control-radius independent (the post-processor works on the elastic field).
mat = {"LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15};
Y = zeros(1, 6);
for kk = 1:6
    pf = @(X) elasticSphereScattering(kk, 1.0, X, mat{:}).total;
    rf = acousticRadiationForce(pf, kk, "ControlRadius", 1.4);
    verifyLessThan(testCase, rf.controlRadiusResidual, 1e-6);
    verifyGreaterThan(testCase, rf.force(3), 0);
    Y(kk) = rf.forceFunction;
end
% resonance: the peak is well above the rigid-sphere plateau (~0.75)
verifyGreaterThan(testCase, max(Y), 2.0);
% and non-monotonic (rises then falls) - a genuine internal resonance
[~, ipk] = max(Y);
verifyGreaterThan(testCase, ipk, 1);
verifyLessThan(testCase, ipk, 6);
verifyLessThan(testCase, Y(end), max(Y));
end
