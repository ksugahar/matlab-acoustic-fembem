function tests = testFsiInvariants
%TESTFSIINVARIANTS Geometry-independent physical invariants of the FSI coupling.
%
% Analytic references (Faran) pin the coupling on the sphere; these invariants
% must hold on ANY geometry/material, so they guard the coupling where no
% analytic solution exists.  Fast lane: acoustic reciprocity of the far-field
% scattering amplitude (two DtN solves).  The heavier optical-theorem energy
% balance (angular far-field integral) lives in validation_test.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testFarFieldScatteringReciprocity(testCase)
% Acoustic reciprocity of the scattering amplitude: f(obs = b, inc = a) equals
% f(obs = -a, inc = -b).  A property of the continuous coupled problem that the
% symmetric FSI block system must reproduce -- independent of any analytic solution.
volFile = fixtureBall();
mat = {"LongitudinalSpeed", 1.6, "ShearSpeed", 0.9, "DensityRatio", 1.15};
k = 1.5;
aHat = [0 0 1];                                   % incidence direction a
bHat = [sin(0.7) 0 cos(0.7)];                     % observation direction b
R = 40;                                           % far-field radius (kR >> 1)

model = FemBemModel(volFile);
fAB = farAmplitude(model, k, mat, aHat, bHat, R);     % inc a, observe b
fBA = farAmplitude(model, k, mat, -bHat, -aHat, R);   % inc -b, observe -a

verifyLessThan(testCase, abs(fAB - fBA) / abs(fAB), 0.05);
end


function f = farAmplitude(model, k, mat, incHat, obsHat, R)
incHat = incHat / norm(incHat);
obsHat = obsHat / norm(obsHat);
inc = struct( ...
    "value", @(X) exp(1i * k * (X * incHat.')), ...
    "grad",  @(X) 1i * k * incHat .* exp(1i * k * (X * incHat.')));
sol = fsiCoupledSolve(model, "Wavenumber", k, mat{:}, ...
    "ExteriorMethod", "dtn", "Incident", inc);
f = R * exp(-1i * k * R) * sol.scatteredAt(R * obsHat);   % p_s ~ f e^{ikR}/R
end


function volFile = fixtureBall()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_ball_maxh018.vol"));
end
