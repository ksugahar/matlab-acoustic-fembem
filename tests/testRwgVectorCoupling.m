function tests = testRwgVectorCoupling
%TESTRWGVECTORCOUPLING H(curl)/RWG vector coupling, ladder stage 6.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testAdjacentTetIndicesAreTheFaceOwners(testCase)
% regression for the ismember-argument bug: every boundary triangle must
% point at the tetrahedron that actually contains its three nodes.
model = FemBemModel(fixturePath("four_tet_interior_node.vol"));
adj = model.mesh.boundaryOrientation.adjacentTetIndices;
verifyEqual(testCase, numel(unique(adj)), size(model.mesh.tet, 1));
for t = 1:size(model.mesh.tri, 1)
    tet = model.mesh.tet(adj(t), :);
    verifyTrue(testCase, all(ismember(model.mesh.tri(t, :), tet)), ...
        sprintf("triangle %d not contained in its adjacent tet", t));
end
end


function testRotatedTangentialTraceIsTheRwgFunction(testCase)
% the FEEC coupling identity, pointwise to machine precision:
%   n x (tangential trace of N_E) = gamma * f_e / l_e on each triangle
% with gamma = -signOut * triEdgeSigns * sigma_pm.
for fixture = ["four_tet_interior_node.vol", "unit_sphere_coarse.vol"]
    model = FemBemModel(fixturePath(fixture));
    maxErr = traceIdentityMaxError(model);
    verifyLessThan(testCase, maxErr, 1e-12, fixture);
end
end


function testRotatedTraceMapReproducesTheField(testCase)
% dof-level contract: for u = sum alpha_E N_E the RWG expansion with
% coefficients C*alpha equals n x u|_Gamma pointwise.
model = FemBemModel(fixturePath("unit_sphere_coarse.vol"));
rng(3);
alpha = randn(model.hcurl.ndof(), 1);
C = model.rwg.rotatedTraceMap(model.hcurl);
c = C * alpha;

quad = SurfaceQuadrature(model.surface, 3);
[Bx, By, Bz] = model.rwg.basisAtQuadrature(quad);
rwgField = [Bx * c, By * c, Bz * c];
nedField = rotatedNedelecTraceAt(model, quad);
traceField = zeros(size(rwgField));
for E = 1:model.hcurl.ndof()
    traceField = traceField + alpha(E) * squeeze(nedField(:, E, :));
end
verifyEqual(testCase, rwgField, traceField, "AbsTol", 1e-11 * max(1, max(abs(alpha))));
end


function testGramIsExactAndPositiveDefinite(testCase)
surface = fixtureSurface("unit_sphere_coarse.vol");
rwg = RwgSpace(surface);
[G3, detail] = rwg.gram();
quad7 = SurfaceQuadrature(surface, 7);
[Bx, By, Bz] = rwg.basisAtQuadrature(quad7);
W = spdiags(quad7.weights, 0, quad7.nPoints(), quad7.nPoints());
G7 = Bx.' * W * Bx + By.' * W * By + Bz.' * W * Bz;
verifyEqual(testCase, full(G3), full(G7), "AbsTol", 1e-12);
verifyEqual(testCase, full(G3), full(G3.'), "AbsTol", 1e-14);
verifyGreaterThan(testCase, min(eig(full(G3))), 0);
verifyEqual(testCase, detail.quadratureOrder, 3);
end


function testVectorSingleLayerIsSymmetricPositiveDefinite(testCase)
surface = fixtureSurface("unit_sphere_coarse.vol");
rwg = RwgSpace(surface);
L = RwgSingleLayer(rwg);
% one-sided semi-analytic assembly leaves test-quadrature asymmetry
% (measured 1.7e-3 at the 3-point rule, same class as the P1 operators)
asym = norm(L.matrix - L.matrix.', "fro") / norm(L.matrix, "fro");
verifyLessThan(testCase, asym, 5e-3);
verifyGreaterThan(testCase, min(eig((L.matrix + L.matrix.') / 2)), 0);
end


function testUniformMagnetizationSphereVectorPotential(testCase)
% K = z_hat x n on the unit sphere is the surface current of a uniformly
% magnetized sphere (M = z_hat, mu0 = 1):
%   A = (1/3) z_hat x x        inside,
%   A = (1/3) z_hat x x / r^3  outside.
% Bands locked from the 2026-07-03 measurement (coarse 3.1%, fine 1.3%;
% geometry-faceting dominated).
cases = struct("fixture", ["unit_sphere_coarse.vol", "unit_sphere_fine.vol"], ...
    "band", [0.05, 0.03]);
for k = 1:2
    surface = fixtureSurface(cases.fixture(k));
    rwg = RwgSpace(surface);
    L = RwgSingleLayer(rwg);
    c = projectRotatedZ(rwg);

    pts = [0 0 0; 0.3 0.2 -0.1; 3 0 0; 0 0 5];
    A = L.vectorPotentialAt(c, pts);
    Aexact = [cross(repmat([0 0 1], 2, 1), pts(1:2, :), 2) / 3
              cross(repmat([0 0 1], 2, 1), pts(3:4, :), 2) ./ (3 * sum(pts(3:4, :).^2, 2).^1.5)];
    relErr = max(abs(A - Aexact), [], "all") / max(abs(Aexact), [], "all");
    verifyLessThan(testCase, relErr, cases.band(k), cases.fixture(k));
end
end


% ---------------------------------------------------------------- helpers

function maxErr = traceIdentityMaxError(model)
s = model.surface;
rwg = model.rwg;
signsOut = s.orientation.triangleOrientationSignsToOutward(:);
adjTet = s.orientation.adjacentTetIndices(:);
areas = s.areas();
quadBary = [4 1 1; 1 4 1; 1 1 4] / 6;
maxErr = 0;
for t = 1:size(s.tri, 1)
    tet = model.mesh.tet(adjTet(t), :);
    X = model.mesh.vtx(tet, :);
    D = [ones(4, 1), X];
    coeff = D \ eye(4);
    gradLambda = coeff(2:4, :).';
    P = s.vtx(s.tri(t, :), :);
    pts = quadBary * P;
    nStored = cross(P(2, :) - P(1, :), P(3, :) - P(1, :));
    nHat = signsOut(t) * nStored / norm(nStored);
    for le = 1:3
        e = rwg.triEdges(t, le);
        d = find(rwg.dofEdgeIds == e, 1);
        if isempty(d), continue, end
        aG = rwg.dofEdgesGlobal(d, 1); bG = rwg.dofEdgesGlobal(d, 2);
        ia = find(tet == aG, 1); ib = find(tet == bG, 1);
        lam = [ones(size(pts, 1), 1), pts] / D;
        N = lam(:, ia) .* gradLambda(ib, :) - lam(:, ib) .* gradLambda(ia, :);
        rotated = cross(repmat(nHat, size(N, 1), 1), N, 2);
        slot = find(rwg.edgeTriangles(e, :) == t, 1);
        sigma = 3 - 2 * slot;
        pOpp = s.vtx(rwg.oppositeVerticesLocal(e, slot), :);
        f = sigma * rwg.dofEdgeLengths(d) / (2 * areas(t)) * (pts - pOpp);
        predicted = -signsOut(t) * rwg.triEdgeSigns(t, le) * sigma ...
            * f / rwg.dofEdgeLengths(d);
        maxErr = max(maxErr, max(abs(rotated - predicted), [], "all"));
    end
end
end


function nedField = rotatedNedelecTraceAt(model, quad)
% n x N_E at the quadrature points, per volume edge: (nPts x nEdges x 3)
s = model.surface;
signsOut = s.orientation.triangleOrientationSignsToOutward(:);
adjTet = s.orientation.adjacentTetIndices(:);
nedField = zeros(quad.nPoints(), model.hcurl.ndof(), 3);
for t = 1:size(s.tri, 1)
    tet = model.mesh.tet(adjTet(t), :);
    X = model.mesh.vtx(tet, :);
    D = [ones(4, 1), X];
    coeff = D \ eye(4);
    gradLambda = coeff(2:4, :).';
    P = s.vtx(s.tri(t, :), :);
    nStored = cross(P(2, :) - P(1, :), P(3, :) - P(1, :));
    nHat = signsOut(t) * nStored / norm(nStored);
    ptsIdx = find(quad.triangleIndex == t);
    pts = quad.points(ptsIdx, :);
    lam = [ones(size(pts, 1), 1), pts] / D;
    localPairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    for k = 1:6
        pair = tet(localPairs(k, :));
        [~, order] = sort(pair);
        ij = localPairs(k, order);              % ascending global orientation
        E = model.hcurl.tetEdges(adjTet(t), k);
        N = lam(:, ij(1)) .* gradLambda(ij(2), :) ...
            - lam(:, ij(2)) .* gradLambda(ij(1), :);
        rotated = cross(repmat(nHat, size(N, 1), 1), N, 2);
        nedField(ptsIdx, E, :) = nedField(ptsIdx, E, :) ...
            + reshape(rotated, [], 1, 3);
    end
end
end


function c = projectRotatedZ(rwg)
% Galerkin projection of K = z_hat x n_out onto the RWG space.
s = rwg.surface;
signsOut = s.orientation.triangleOrientationSignsToOutward(:);
quad = SurfaceQuadrature(s, 3);
[Bx, By, ~] = rwg.basisAtQuadrature(quad);
W = spdiags(quad.weights, 0, quad.nPoints(), quad.nPoints());
nrmStored = cross(s.vtx(s.tri(:, 2), :) - s.vtx(s.tri(:, 1), :), ...
    s.vtx(s.tri(:, 3), :) - s.vtx(s.tri(:, 1), :), 2);
nrmOut = signsOut .* nrmStored ./ sqrt(sum(nrmStored.^2, 2));
nAtPts = nrmOut(quad.triangleIndex, :);
b = Bx.' * (W * (-nAtPts(:, 2))) + By.' * (W * nAtPts(:, 1));
[G, ~] = rwg.gram();
c = G \ b;
end


function surface = fixtureSurface(name)
mesh = VolMesh(fixturePath(name));
surface = mesh.boundary();
end


function path = fixturePath(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
path = string(fullfile(repoRoot, "fixtures", "mesh_topology", name));
end
