function fem = assembleNedelec0TetMatrices(model, materialCoef)
%ASSEMBLENEDELEC0TETMATRICES Assemble readable first-order HCurl matrices.
%
% The local basis is N_ij = lambda_i grad(lambda_j) -
% lambda_j grad(lambda_i), one basis function per oriented tetrahedron edge.

arguments
    model (1,1) struct
    materialCoef double = 1
end

topology = buildFirstOrderTopology(model);
nodes = model.lukas.geo.nodes;
tet = model.lukas.geo.conn_matrix;
nEdges = size(topology.hcurl.edges, 1);
nTets = size(tet, 1);
localPairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

if isscalar(materialCoef)
    materialCoef = repmat(materialCoef, nTets, 1);
else
    materialCoef = materialCoef(:);
end
if numel(materialCoef) ~= nTets
    error("assembleNedelec0TetMatrices:material", "materialCoef must be scalar or one value per tetrahedron.");
end

ii = zeros(36 * nTets, 1);
jj = zeros(36 * nTets, 1);
massVals = zeros(36 * nTets, 1);
curlVals = zeros(36 * nTets, 1);
localMass = zeros(6, 6, nTets);
localCurlCurl = zeros(6, 6, nTets);
localCurls = zeros(6, 3, nTets);
volumes = zeros(nTets, 1);

cursor = 1;
for e = 1:nTets
    X = nodes(tet(e, :), :);
    [Me, Ce, curls, volume] = localNedelec0(X, materialCoef(e), localPairs);
    signs = topology.hcurl.tetEdgeSigns(e, :).';
    Me = (signs * signs.') .* Me;
    Ce = (signs * signs.') .* Ce;
    curls = signs .* curls;

    localMass(:, :, e) = Me;
    localCurlCurl(:, :, e) = Ce;
    localCurls(:, :, e) = curls;
    volumes(e) = volume;

    ids = topology.hcurl.tetEdges(e, :);
    [I, J] = ndgrid(ids, ids);
    span = cursor:(cursor + 35);
    ii(span) = I(:);
    jj(span) = J(:);
    massVals(span) = Me(:);
    curlVals(span) = Ce(:);
    cursor = cursor + 36;
end

fem = struct();
fem.family = "Nedelec0";
fem.cell = "tetrahedron";
fem.edges = topology.hcurl.edges;
fem.mass = sparse(ii, jj, massVals, nEdges, nEdges);
fem.curlCurl = sparse(ii, jj, curlVals, nEdges, nEdges);
fem.localMass = localMass;
fem.localCurlCurl = localCurlCurl;
fem.localCurls = localCurls;
fem.volumes = volumes;
fem.materialCoef = materialCoef;
end


function [Me, Ce, curls, volume] = localNedelec0(X, materialCoef, localPairs)
D = [ones(4, 1), X];
detD = det(D);
volume = abs(detD) / 6;
if volume <= eps
    error("assembleNedelec0TetMatrices:degenerate", "Degenerate tetrahedron with near-zero volume.");
end

coeff = D \ eye(4);
gradLambda = coeff(2:4, :).';
lambdaMass = volume / 20 * (ones(4, 4) + eye(4));

Me = zeros(6, 6);
Ce = zeros(6, 6);
curls = zeros(6, 3);
for a = 1:6
    i = localPairs(a, 1);
    j = localPairs(a, 2);
    curls(a, :) = 2 * cross(gradLambda(i, :), gradLambda(j, :));
    for b = 1:6
        k = localPairs(b, 1);
        l = localPairs(b, 2);
        Me(a, b) = materialCoef * ( ...
            dot(gradLambda(j, :), gradLambda(l, :)) * lambdaMass(i, k) ...
            - dot(gradLambda(j, :), gradLambda(k, :)) * lambdaMass(i, l) ...
            - dot(gradLambda(i, :), gradLambda(l, :)) * lambdaMass(j, k) ...
            + dot(gradLambda(i, :), gradLambda(k, :)) * lambdaMass(j, l));
        Ce(a, b) = materialCoef * volume * dot(curls(a, :), ...
            2 * cross(gradLambda(k, :), gradLambda(l, :)));
    end
end
end
