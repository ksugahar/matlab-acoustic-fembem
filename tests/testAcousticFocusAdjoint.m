function tests = testAcousticFocusAdjoint
%TESTACOUSTICFOCUSADJOINT Adjoint AD through the rigid BEM solve (wavefront).
%
% Locks reverse-mode automatic differentiation for the phased-array
% focusing problem: the gradient of the focused intensity J = |u(target)|^2
% with respect to the complex source amplitudes, computed by ONE adjoint
% (transpose) solve through the rigid-scattering BEM system, matches central
% finite differences (2026-07-04 measurements: forward affine residual
% 4.4e-18, adjoint-vs-FD 1.7e-10). A short gradient ascent must focus energy
% at the target (the wavefront-synthesis proof). This is the seed of the
% acoustic radiation-force / thrust design lane.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function [surface, sources, target, k] = focusProblem()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
mesh = VolMesh(fullfile(repoRoot, "fixtures", "mesh_topology", ...
    "unit_sphere_fine.vol"));
surface = mesh.boundary();
k = 2.0;
n = 8;
ang = (0:n-1).' / n * 2*pi;
sources = [2.5*cos(ang), 2.5*sin(ang), -3.0*ones(n, 1)];
target = [0 0 2.5];   % behind the sphere from the array
end


function testForwardMapIsExactlyAffine(testCase)
% u(target) must equal w * amplitudes to machine precision - this validates
% the adjoint row w independently of any finite difference.
[surface, sources, target, k] = focusProblem();
rng(3);
p = randn(8, 1) + 1i*randn(8, 1);
res = acousticFocusAdjoint(surface, sources, target, k, p);
verifyLessThan(testCase, res.forwardLinearityResidual, 1e-10);
verifyEqual(testCase, res.field, res.sensitivityRow * p, "AbsTol", 1e-10);
verifyEqual(testCase, res.adjointSolves, 1);
end


function testAdjointGradientMatchesFiniteDifference(testCase)
% the reverse-mode payoff: one adjoint solve for the whole gradient, equal
% to central finite differences (which cost one forward re-solve each).
[surface, sources, target, k] = focusProblem();
rng(7);
p = randn(8, 1) + 1i*randn(8, 1);
res = acousticFocusAdjoint(surface, sources, target, k, p, ...
    "GradientCheck", true);
verifyEqual(testCase, res.status, "ok");
verifyTrue(testCase, res.checks.gradientMatchesFiniteDifference);
verifyLessThan(testCase, res.gradientCheckRelError, 1e-6);   % measured 1.7e-10
end


function testAscentDirectionIsCorrectAndFocuses(testCase)
% wavefront synthesis. The forward map u = w * p is fixed (w = sensitivity
% row), so the ascent is a pure-arithmetic loop with NO further BEM solves -
% one adjoint call gives w, the rest is free. Rigorous, seed-independent
% correctness: a small step along +ascentDirection RAISES J and -ascent
% LOWERS it (the sign that the wrong Wirtinger derivative, dJ/dp instead of
% dJ/dconj(p), would flip - a bug near-zero fields hide but this start,
% away from the minimum, exposes).
[surface, sources, target, k] = focusProblem();
rng(11);
p = randn(8, 1) + 1i*randn(8, 1);
res = acousticFocusAdjoint(surface, sources, target, k, p);
w = res.sensitivityRow;

Jof = @(pp) abs(w * pp)^2;
ascentAt = @(pp) ascentDir(w, pp);

d = ascentAt(p);
alpha = 1e-3 / max(abs(d));
verifyGreaterThan(testCase, Jof(p + alpha * d), res.objective);   % +ascent up
verifyLessThan(testCase, Jof(p - alpha * d), res.objective);      % -ascent down

% a damped ascent must monotonically focus energy (measured factor >> 1)
J0 = Jof(p);
Jprev = J0;
monotone = true;
for it = 1:200
    d = ascentAt(p);
    p = p + (5e-3 / max(abs(d))) * d;
    Jc = Jof(p);
    if Jc < Jprev - 1e-12
        monotone = false;
    end
    Jprev = Jc;
end
verifyTrue(testCase, monotone);
verifyGreaterThan(testCase, Jc / J0, 5);
end


function d = ascentDir(w, p)
% steepest-ascent direction of J = |w p|^2: dJ/dconj(p) = 2 u conj(w).
u = w * p;
d = (2 * u * conj(w)).';
end
