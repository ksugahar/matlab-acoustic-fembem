function tests = testLaplacePanelIntegrals
%TESTLAPLACEPANELINTEGRALS Analytic 1/r and P1/r triangle integrals.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testMatchesSubdivisionQuadratureOffPlane(testCase)
rng(7);
for t = 1:5
    V = randn(3, 3);
    while norm(cross(V(2,:)-V(1,:), V(3,:)-V(1,:))) < 0.3
        V = randn(3, 3);
    end
    ctr = mean(V, 1);
    nrm = cross(V(2,:)-V(1,:), V(3,:)-V(1,:)); nrm = nrm / norm(nrm);
    diam = max(vecnorm(V - circshift(V, 1), 2, 2));
    pts = [ctr + 2.5*diam*nrm; ctr + 0.1*diam*nrm; V(1,:) + 0.4*diam*nrm];

    [I0a, I1a] = laplacePanelIntegrals(V, pts);
    [I0n, I1n] = referenceBySubdivision(V, pts, 6);
    verifyEqual(testCase, I0a, I0n, "RelTol", 1e-7);
    verifyEqual(testCase, I1a, I1n, "AbsTol", 1e-7 * max(abs(I0n)));
end
end


function testSelfPointMatchesPolarIntegration(testCase)
% observation point inside the triangle, in plane: the singular case the
% Galerkin operator hits on every diagonal panel.
V = [0 0 0; 1.3 0 0; 0.4 0.9 0];
bary = [1/3 1/3 1/3; 0.6 0.3 0.1];
pts = bary * V;

[I0, I1] = laplacePanelIntegrals(V, pts);

% polar reference for the constant integral: split into 3 apex triangles
for q = 1:size(pts, 1)
    x0 = pts(q, :);
    I0ref = 0;
    for e = 1:3
        a = V(e, :); b = V(mod(e, 3) + 1, :);
        u = a - x0; v = b - x0;
        phiMax = atan2(norm(cross(u, v)), dot(u, v));
        uh = u / norm(u);
        wh = cross(cross(u, v), u); wh = wh / norm(wh);
        rayLen = @(phi) arrayfun(@(ph) ...
            edgeHit(x0, cos(ph)*uh + sin(ph)*wh, u, v), phi);
        I0ref = I0ref + integral(rayLen, 0, phiMax, "AbsTol", 1e-13);
    end
    verifyEqual(testCase, I0(q), I0ref, "RelTol", 1e-10);
end
verifyEqual(testCase, sum(I1, 2), I0, "AbsTol", 1e-12 * max(abs(I0)));
end


function testPartitionOfUnity(testCase)
rng(11);
V = randn(3, 3) + [3 0 0; 0 0 0; 0 0 0];
pts = [randn(4, 3); mean(V, 1)];
[I0, I1] = laplacePanelIntegrals(V, pts);
verifyEqual(testCase, sum(I1, 2), I0, "AbsTol", 1e-12 * max(abs(I0)));
end


function testRejectsDegenerateTriangle(testCase)
V = [0 0 0; 1 0 0; 2 0 0];
verifyError(testCase, @() laplacePanelIntegrals(V, [0 1 0]), ...
    "laplacePanelIntegrals:degenerate");
end


function L = edgeHit(x0, dirVec, u, v)
A = [dirVec.', (u - v).'];
sol = A \ u.';
L = sol(1);
end


function [I0, I1] = referenceBySubdivision(V, pts, levels)
tris = {V};
for L = 1:levels
    next = cell(1, 4 * numel(tris));
    for i = 1:numel(tris)
        T = tris{i};
        m12 = (T(1,:)+T(2,:))/2; m23 = (T(2,:)+T(3,:))/2; m31 = (T(3,:)+T(1,:))/2;
        next{4*i-3} = [T(1,:); m12; m31];
        next{4*i-2} = [m12; T(2,:); m23];
        next{4*i-1} = [m31; m23; T(3,:)];
        next{4*i}   = [m12; m23; m31];
    end
    tris = next;
end
a1 = (6 - sqrt(15)) / 21;  a2 = (6 + sqrt(15)) / 21;
w1 = (155 - sqrt(15)) / 1200; w2 = (155 + sqrt(15)) / 1200;
bq = [1/3 1/3 1/3
      a1 a1 1-2*a1; a1 1-2*a1 a1; 1-2*a1 a1 a1
      a2 a2 1-2*a2; a2 1-2*a2 a2; 1-2*a2 a2 a2];
wq = [9/40; w1; w1; w1; w2; w2; w2];

n = cross(V(2,:)-V(1,:), V(3,:)-V(1,:));
A2 = norm(n); nh = n / A2;
idx = [2 3; 3 1; 1 2];

I0 = zeros(size(pts, 1), 1);
I1 = zeros(size(pts, 1), 3);
for i = 1:numel(tris)
    T = tris{i};
    area = 0.5 * norm(cross(T(2,:)-T(1,:), T(3,:)-T(1,:)));
    gp = bq * T;
    lam = zeros(size(gp, 1), 3);
    for k = 1:3
        a = V(idx(k, 1), :); b = V(idx(k, 2), :);
        for r = 1:size(gp, 1)
            lam(r, k) = dot(cross(b - a, gp(r, :) - a), nh) / A2;
        end
    end
    for p = 1:size(pts, 1)
        r = sqrt(sum((gp - pts(p, :)).^2, 2));
        I0(p) = I0(p) + area * sum(wq ./ r);
        I1(p, :) = I1(p, :) + area * (wq ./ r).' * lam;
    end
end
end
