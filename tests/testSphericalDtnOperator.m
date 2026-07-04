function tests = testSphericalDtnOperator
%TESTSPHERICALDTNOPERATOR Exact spherical Helmholtz DtN (Kelvin on the sphere).
%
% The low-rank exterior FAST PATH that replaces the dense Galerkin single/
% double layer when the acoustic truncation surface is a sphere - the operator
% the Kelvin transformation (and its radiating extension, Sugahara IEICE 2024)
% represents on the sphere. Gated from the 2026-07-04 measurements:
%   sphere fit deviation 2e-16, Gram cond 7.4 (well-posed);
%   INDEPENDENT point-source DtN check 2.6e-5 at degree 10 - the operator maps
%     an outgoing field's Dirichlet trace to its ANALYTIC Neumann trace (not the
%     harmonics it is built from), and converges exponentially in the degree
%     (7e-3, 1e-3, 2e-4, 3e-5 at degree 4/6/8/10) to the ~1e-5 mesh/P1 floor;
%   fail-loud on a non-spherical surface (no silent fallback to the sphere).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function surface = fixtureSurface(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
surface = VolMesh(fullfile(repoRoot, "fixtures", "mesh_topology", name)).boundary();
end


function [pG, qG] = pointSourceTraces(surface, k, xs, center)
% Dirichlet + analytic Neumann (outward) traces of an interior point source's
% outgoing exterior field on the sphere of the given center, for the
% independent DtN check (the outward normal is measured from that center).
vtx = surface.vtx;
dirs = (vtx - center) ./ sqrt(sum((vtx - center).^2, 2));
rv = vtx - xs; rr = sqrt(sum(rv.^2, 2));
pG = exp(1i * k * rr) ./ (4 * pi * rr);
qG = pG .* (1i * k - 1 ./ rr) .* (sum(rv .* dirs, 2) ./ rr);
end


function testSphereGateAndHealth(testCase)
op = sphericalDtnOperator(fixtureSurface("unit_sphere_fine.vol"), "Wavenumber", 2.0);
verifyEqual(testCase, op.status, "ok");
verifyLessThan(testCase, op.sphericityDeviation, 1e-6);   % nodes on the unit sphere
verifyLessThan(testCase, op.gramCondition, 100);          % measured 7.4
verifyEqual(testCase, op.radius, 1.0, "AbsTol", 1e-6);
end


function testIndependentPointSourceDtN(testCase)
% an off-center INTERIOR point source radiates an outgoing exterior field; the
% DtN must map its Dirichlet trace to its analytic Neumann trace on the sphere.
k = 2.0;
surface = fixtureSurface("unit_sphere_fine.vol");
op = sphericalDtnOperator(surface, "Wavenumber", k);
[pG, qG] = pointSourceTraces(surface, k, [0.3 0.1 -0.2], op.center);
errA = norm(op.apply(pG) - qG) / norm(qG);
verifyLessThan(testCase, errA, 1e-4);                     % measured 2.6e-5 (degree 10)
end


function testDegreeConverged(testCase)
% the DtN is exact per multipole, so the point-source check CONVERGES in the
% degree (measured 7e-3, 1e-3, 2e-4, 3e-5 at degree 4/6/8/10, down to the
% ~1e-5 mesh/P1 floor) - higher degree is strictly better until the floor.
k = 2.0;
surface = fixtureSurface("unit_sphere_fine.vol");
op6 = sphericalDtnOperator(surface, "Wavenumber", k, "Degree", 6);
op10 = sphericalDtnOperator(surface, "Wavenumber", k, "Degree", 10);
[pG, qG] = pointSourceTraces(surface, k, [0.3 0.1 -0.2], op10.center);
e6 = norm(op6.apply(pG) - qG) / norm(qG);
e10 = norm(op10.apply(pG) - qG) / norm(qG);
verifyLessThan(testCase, e6, 5e-3);                       % measured 1.1e-3
verifyLessThan(testCase, e10, e6);                        % higher degree strictly better
verifyLessThan(testCase, e10, 1e-4);                      % measured 2.6e-5
end


function testFailLoudOnNonSphere(testCase)
% a genuinely non-spherical surface (a 5-sphere chain) must RAISE, not silently
% fall back - the DtN fast path is only valid on a sphere truncation.
verifyError(testCase, ...
    @() sphericalDtnOperator(fixtureSurface("soft_sphere_chain_5.vol"), "Wavenumber", 1.0), ...
    "sphericalDtnOperator:notSpherical");
end
