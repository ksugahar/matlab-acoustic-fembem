classdef Nedelec0Space
%NEDELEC0SPACE Lowest-order HCurl (Nedelec0) edge space on tetrahedra.
%
%   space = Nedelec0Space(mesh);
%   [M, detail] = space.mass(materialCoef);       % int c N_a . N_b dx
%   [C, detail] = space.curlCurl(materialCoef);   % int c curl(N_a) . curl(N_b) dx
%
% One dof per oriented volume edge. The local basis on each tetrahedron is
%
%   N_ij = lambda_i grad(lambda_j) - lambda_j grad(lambda_i)
%
% with the global orientation fixed by sorting node ids: an edge stored as
% (i, j) with i < j is positive, and tetEdgeSigns records the per-tet flip.

properties (Constant)
    family = "Nedelec"
    order = 0
    cell = "tetrahedron"
    basis = "Nedelec0"
    localEdgePairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
end

properties
    mesh          % VolMesh carrying vtx/tet
    edges         % oriented volume edges (nEdges x 2), sorted node ids
    tetEdges      % edge id per tetrahedron local edge (nTets x 6)
    tetEdgeSigns  % +1/-1 local-vs-global orientation (nTets x 6)
end

methods
    function space = Nedelec0Space(mesh)
        arguments
            mesh (1,1) VolMesh
        end
        space.mesh = mesh;
        [space.edges, space.tetEdges, space.tetEdgeSigns] = ...
            tetEdgeTopology(mesh.tet, space.localEdgePairs);
    end

    function n = ndof(space)
        n = size(space.edges, 1);
    end

    function [M, detail] = mass(space, materialCoef)
        %MASS Assemble int c N_a . N_b dx over the oriented edge dofs.
        arguments
            space (1,1) Nedelec0Space
            materialCoef double = 1
        end
        [M, ~, detail] = assembleEdgeMatrices(space, materialCoef);
    end

    function [C, detail] = curlCurl(space, materialCoef)
        %CURLCURL Assemble int c curl(N_a) . curl(N_b) dx.
        arguments
            space (1,1) Nedelec0Space
            materialCoef double = 1
        end
        [~, C, detail] = assembleEdgeMatrices(space, materialCoef);
    end

    function [M, C, detail] = matrices(space, materialCoef)
        %MATRICES Assemble mass and curl-curl in one per-tet sweep.
        arguments
            space (1,1) Nedelec0Space
            materialCoef double = 1
        end
        [M, C, detail] = assembleEdgeMatrices(space, materialCoef);
    end
end

methods (Access = private)
    function [M, C, detail] = assembleEdgeMatrices(space, materialCoef)
        nodes = space.mesh.vtx;
        tet = space.mesh.tet;
        nEdges = size(space.edges, 1);
        nTets = size(tet, 1);
        if isscalar(materialCoef)
            materialCoef = repmat(materialCoef, nTets, 1);
        else
            materialCoef = materialCoef(:);
        end
        if numel(materialCoef) ~= nTets
            error("Nedelec0Space:material", ...
                "materialCoef must be scalar or one value per tetrahedron.");
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
            [Me, Ce, curls, volume] = ...
                localNedelec0(X, materialCoef(e), space.localEdgePairs);
            signs = space.tetEdgeSigns(e, :).';
            Me = (signs * signs.') .* Me;
            Ce = (signs * signs.') .* Ce;
            curls = signs .* curls;

            localMass(:, :, e) = Me;
            localCurlCurl(:, :, e) = Ce;
            localCurls(:, :, e) = curls;
            volumes(e) = volume;

            ids = space.tetEdges(e, :);
            [I, J] = ndgrid(ids, ids);
            span = cursor:(cursor + 35);
            ii(span) = I(:);
            jj(span) = J(:);
            massVals(span) = Me(:);
            curlVals(span) = Ce(:);
            cursor = cursor + 36;
        end

        M = sparse(ii, jj, massVals, nEdges, nEdges);
        C = sparse(ii, jj, curlVals, nEdges, nEdges);
        detail = struct( ...
            "family", space.basis, ...
            "cell", space.cell, ...
            "edges", space.edges, ...
            "mass", M, ...
            "curlCurl", C, ...
            "localMass", localMass, ...
            "localCurlCurl", localCurlCurl, ...
            "localCurls", localCurls, ...
            "volumes", volumes, ...
            "materialCoef", materialCoef);
    end
end
end


function [edges, tetEdges, tetEdgeSigns] = tetEdgeTopology(tet, localPairs)
%TETEDGETOPOLOGY Unique oriented volume edges with per-tet signs.

nTets = size(tet, 1);
raw = zeros(6 * nTets, 2);
signsRaw = zeros(6 * nTets, 1);

for e = 1:nTets
    for k = 1:6
        row = (e - 1) * 6 + k;
        pair = tet(e, localPairs(k, :));
        sortedPair = sort(pair);
        raw(row, :) = sortedPair;
        if isequal(pair, sortedPair)
            signsRaw(row) = 1;
        else
            signsRaw(row) = -1;
        end
    end
end

[edges, ~, ic] = unique(raw, "rows");
tetEdges = reshape(ic, 6, nTets).';
tetEdgeSigns = reshape(signsRaw, 6, nTets).';
end


function [Me, Ce, curls, volume] = localNedelec0(X, materialCoef, localPairs)
%LOCALNEDELEC0 Per-tet 6x6 mass and curl-curl for N_ij basis functions.

D = [ones(4, 1), X];
detD = det(D);
volume = abs(detD) / 6;
if volume <= eps
    error("Nedelec0Space:degenerate", "Degenerate tetrahedron with near-zero volume.");
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
