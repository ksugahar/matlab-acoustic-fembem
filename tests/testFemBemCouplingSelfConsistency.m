function tests = testFemBemCouplingSelfConsistency
%TESTFEMBEMCOUPLINGSELFCONSISTENCY Structure-FEM / acoustic-BEM coupling self-checks.
%
% The FEM-BEM coupling verifies by a MODULAR self-consistency decomposition, no
% analytic reference:
%
%   (a) acoustic BEM operators  -> testBemLayerSelfConsistency (Calderon/Green)
%   (b) elastic FEM operator    -> the 6 rigid-body modes are strain-free
%                                  (null space of K).  Works for BOTH the
%                                  VOLUME-element (elasticityMatrices) and the
%                                  SURFACE-element/shell (surfaceElasticityMatrices)
%                                  variant -- the same patch test, swapped operator.
%   (c) interface geometry      -> the closed wetted surface has zero total
%                                  vector area, so a rigid translation injects no
%                                  net monopole into the fluid.
%   (d) assembled coupled solve -> the scattered boundary traces (p_s, q_s) from
%                                  fsiCoupledSolve are consistent exterior
%                                  radiating Cauchy data: D[p_s] - S[q_s] VANISHES
%                                  inside the scatterer (the null-field).
%
% (b)+(c) are variant-agnostic in structure (rigidBodyModes below is shared), so
% a surface-element shell FEM plugs into the identical coupling harness; only the
% elastic operator and -- for a two-sided-wetted shell -- the exterior BEM (open
% surface + hypersingular W) change.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function coarse = coarseVol()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
coarse = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_coarse.vol");
end


function U = rigidBodyModes(vtx)
%RIGIDBODYMODES The 6 rigid-body displacement modes (3 translations + 3
% rotations) on a node set, in the (node n, component c) at 3*(n-1)+c layout
% shared by elasticityMatrices, surfaceElasticityMatrices, and the coupling G.
nV = size(vtx, 1);
U = zeros(3*nV, 6);
U(1:3:end, 1) = 1; U(2:3:end, 2) = 1; U(3:3:end, 3) = 1;      % translations
for a = 1:3
    e = zeros(1, 3); e(a) = 1;
    rot = cross(repmat(e, nV, 1), vtx, 2);                    % omega x r
    U(1:3:end, 3+a) = rot(:, 1);
    U(2:3:end, 3+a) = rot(:, 2);
    U(3:3:end, 3+a) = rot(:, 3);
end
end


function r = rigidNullResidual(K, U)
%RIGIDNULLRESIDUAL max_i ||K u_i|| / (||K|| ||u_i||) over the rigid modes U.
Kn = normest(K);
r = 0;
for i = 1:size(U, 2)
    r = max(r, norm(K*U(:, i)) / (Kn * norm(U(:, i))));
end
end


function testVolumeElasticStiffnessAnnihilatesRigidModes(testCase)
% (b) VOLUME-element variant: vector P1 tet elasticity K annihilates all 6
% rigid-body modes; the consistent mass M stores positive kinetic energy.
model = FemBemModel(coarseVol());
vtx = model.mesh.vtx;
[K, M] = elasticityMatrices(model.mesh, 1.0, 0.7, 1.2);
U = rigidBodyModes(vtx);
verifyLessThan(testCase, rigidNullResidual(K, U), 1e-10);     % measured ~5e-17
verifyGreaterThan(testCase, U(:,1).' * M * U(:,1), 0);        % translation KE > 0
end


function testSurfaceMembraneStiffnessAnnihilatesRigidModes(testCase)
% (b) SURFACE-element variant: the P1 membrane K (surfaceElasticityMatrices)
% annihilates the same 6 rigid modes, while a non-rigid stretch stores positive
% energy -- so the surface-element shell FEM self-verifies by the IDENTICAL
% patch test as the volume variant.
surface = VolMesh(coarseVol()).boundary();
vtx = surface.vtx; nV = size(vtx, 1);
[K, M] = surfaceElasticityMatrices(surface, 1.0, 0.3, 0.05, 1.0);
U = rigidBodyModes(vtx);
verifyLessThan(testCase, rigidNullResidual(K, U), 1e-10);     % measured ~7e-17

us = zeros(3*nV, 1); us(1:3:end) = vtx(:, 1);                 % global x-stretch
verifyGreaterThan(testCase, us.' * K * us, 1e-6 * normest(K) * (us.'*us));
verifyGreaterThan(testCase, U(:,1).' * M * U(:,1), 0);
end


function testClosedWettedSurfaceHasZeroVectorArea(testCase)
% (c) interface geometry: for a CLOSED wetted surface int_Gamma n dS = 0, so a
% rigid translation of the structure injects zero net monopole into the fluid.
surface = VolMesh(coarseVol()).boundary();
signs = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri; vtx = surface.vtx;
vecArea = [0 0 0]; totArea = 0;
for t = 1:size(tri, 1)
    X = vtx(tri(t, :), :);
    cr = cross(X(2, :) - X(1, :), X(3, :) - X(1, :));
    vecArea = vecArea + 0.5 * signs(t) * cr;
    totArea = totArea + 0.5 * norm(cr);
end
verifyLessThan(testCase, norm(vecArea) / totArea, 1e-10);     % measured ~6e-17
end


function testCoupledSolveScatteredTracesAreExteriorRadiating(testCase)
% (d) assembled coupling: solve the elastic-FEM / acoustic-BEM FSI, then confirm
% the scattered boundary traces are consistent exterior radiating Cauchy data --
% D[p_s] - S[q_s] vanishes INSIDE the scatterer (the exterior representation's
% null-field), a self-consistency of the full coupled assembly.
model = FemBemModel(coarseVol());
sol = fsiCoupledSolve(model, "Wavenumber", 2.0, "LongitudinalSpeed", 2.0, ...
    "ShearSpeed", 1.0, "DensityRatio", 1.5, "QuadratureOrder", 3);
verifyEqual(testCase, sol.status, "ok");

interiorPts = [0 0 0.3; 0.2 0 0; 0 0.25 0.1];
exteriorPts = [0 0 3; 2.5 0 0; 0 3 0];
nullField = norm(sol.scatteredAt(interiorPts)) / norm(sol.scatteredAt(exteriorPts));
verifyLessThan(testCase, nullField, 5e-2);                    % measured 2.7e-2
verifyFalse(testCase, isreal(sol.scatteredAt(exteriorPts))); % genuinely radiating
end
