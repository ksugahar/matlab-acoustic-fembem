function tests = testElasticThrustAdjoint
%TESTELASTICTHRUSTADJOINT Elastic-bead wavefront-synthesis thrust adjoint.
%
% The acoustic RADIATION FORCE on a solid ELASTIC bead from a phased array, and
% its Wirtinger gradient, computed THROUGH the FSI coupled solve - the last mile
% of the ultrasonic-thrust story (design the array phases to steer the force on
% an elastic bead, whose internal resonances a rigid/soft bead lacks). Locked
% from the 2026-07-04 measurements:
%   the force quadratic form p^H Q_i p reproduces the vectorised direct Brillouin
%     integral to ~1e-15 (independent assembly - catches an outer-product
%     conjugation error), and that matches the golden acousticRadiationForce;
%   the force is control-radius independent to ~1e-6 (div T = 0);
%   the Wirtinger gradient matches central finite differences to ~1e-13;
%   forceForm is reusable (no re-solve): a step along the ascent raises F_z.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function model = fixtureModel(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
model = FemBemModel(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end


function testThrustForceGradientAndAscent(testCase)
k = 2.0;
model = fixtureModel("unit_ball_maxh018.vol");
nSrc = 4;
ang = (0:nSrc-1).' / nSrc * 2*pi;
sources = [2.5*cos(ang), 2.5*sin(ang), -3*ones(nSrc,1)];   % ring array below the bead
amps = exp(1i * ang * 0.7);

res = elasticThrustAdjoint(model, sources, k, amps, ...
    "NMu", 8, "NPhi", 12, "GradientCheck", true);

verifyEqual(testCase, res.status, "ok");
verifyLessThan(testCase, res.consistencyError, 1e-8);        % Q == direct integral
verifyLessThan(testCase, res.independentForceError, 1e-8);   % direct == golden arf
verifyLessThan(testCase, res.controlRadiusResidual, 5e-3);   % div T = 0 (quadrature-limited at 8/12; 3e-6 at 12/24)
verifyLessThan(testCase, res.gradientCheckRelError, 1e-6);   % Wirtinger vs FD
verifyEqual(testCase, res.fsiSolves, nSrc);

% forceForm reuse (no re-solve): a small step along the ascent raises F_z
Q3 = res.forceForm{3};
d = res.ascentDirection / norm(res.ascentDirection);
Fz0 = real(amps' * Q3 * amps);
Fz1 = real((amps + 1e-3 * d)' * Q3 * (amps + 1e-3 * d));
verifyGreaterThan(testCase, Fz1, Fz0);
end
