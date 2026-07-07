function tests = testElasticBoxScattering
%TESTELASTICBOXSCATTERING Elastic FSI on a rectangular (box) scatterer.
%
% The frequency-domain elastic FSI (fsiCoupledSolve) is validated against the
% analytic Faran series on a SPHERE (testFsiCoupledSolve).  This exercises the
% SAME solver on a non-separable RECTANGULAR box -- the BEM exterior handles
% arbitrary geometry, so no analytic reference exists -- and locks two
% shape-independent physics invariants instead: acoustic reciprocity of the
% far-field scattering amplitude f(x_hat; d) = f(-d; -x_hat), and the Sommerfeld
% 1/r radiation decay of the scattered field.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
volFile = string(fullfile(tempdir, ...
    "elastic_box_scatterer_" + char(java.util.UUID.randomUUID()) + ".vol"));
structuredBoxVol(volFile, Size=[1.5 1.0 0.8], Cells=[6 4 4], ...
    MaterialName="solid", BoundaryName="interface");
testCase.TestData.volFile = volFile;
testCase.TestData.model = FemBemModel(volFile);
end


function teardownOnce(testCase)
if isfile(testCase.TestData.volFile)
    delete(testCase.TestData.volFile);
end
end


function opts = materialOptions()
opts = {"LongitudinalSpeed", 2.0, "ShearSpeed", 1.0, "DensityRatio", 1.5, ...
    "QuadratureOrder", 3};
end


function testBoxFsiSolvesCleanlyAndRadiates(testCase)
% The elastic FSI solves cleanly on a genuine rectangular box and the scattered
% field decays like the Sommerfeld 1/r far field.
model = testCase.TestData.model;
mat = materialOptions;
sol = fsiCoupledSolve(model, "Wavenumber", 2.0, "ExteriorMethod", "bem", ...
    mat{:});
verifyEqual(testCase, sol.status, "ok");
verifyLessThan(testCase, sol.residualNorm, 1e-8);
verifyTrue(testCase, ~isreal(sol.surfacePressure));   % genuine complex scattered field

xd = [1 0 0];
p_r  = sol.scatteredAt(30 * xd);
p_2r = sol.scatteredAt(60 * xd);
verifyEqual(testCase, abs(p_2r) / abs(p_r), 0.5, "AbsTol", 0.02);   % 1/r decay
end


function testBoxScatteringReciprocity(testCase)
% Acoustic far-field reciprocity f_inf(x_hat; d) = f_inf(-d; -x_hat) holds for
% ANY scatterer shape; here it is satisfied to discretisation error (~1.5%).
model = testCase.TestData.model;
k = 2.0;
r = 60;
mat = materialOptions;
planeWave = @(d) struct( ...
    "value", @(X) exp(1i*k*(X*d(:))), ...
    "grad",  @(X) 1i*k*exp(1i*k*(X*d(:))) .* repmat(d(:).', size(X,1), 1));
farAmp = @(sol, xhat) r * exp(-1i*k*r) * sol.scatteredAt(r * xhat(:).');
solveFar = @(d, xhat) farAmp( ...
    fsiCoupledSolve(model, "Wavenumber", k, "ExteriorMethod", "bem", ...
        "Incident", planeWave(d), mat{:}), xhat);

pairs = {[0 0 1], [1 0 0]; [1 1 0]/sqrt(2), [0 0 1]};   % {incident d, observe x_hat}
for p = 1:size(pairs, 1)
    d = pairs{p, 1};
    xhat = pairs{p, 2};
    fForward = solveFar(d, xhat);
    fReverse = solveFar(-xhat, -d);
    verifyLessThan(testCase, abs(fForward - fReverse) / abs(fForward), 3e-2, ...
        sprintf("far-field reciprocity broken for pair %d", p));
end
end
