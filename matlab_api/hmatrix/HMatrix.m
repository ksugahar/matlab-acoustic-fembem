classdef HMatrix
%HMATRIX Readable H-matrix for the Laplace single-layer teaching kernel.
%
%   H = HMatrix(target);                  % points, SurfaceMesh, or FemBemModel
%   H = HMatrix(target, source, "LeafSize", 8, "Eta", 2.0);
%   y = H.matvec(x);                      % or H * x
%   s = H.stats();                        % block-tree compression statistics
%
% This is intentionally small and inspectable. It demonstrates the core
% NGSolve.BEM/Gypsilab ideas without production engineering:
%
%   - cluster trees by recursive bisection of the longest bbox axis
%   - admissible far blocks: diameter <= Eta * distance
%   - low-rank far-field blocks by truncated SVD
%   - dense near-field blocks
%   - a recursive block-tree matvec
%
% The kernel is the collocation single layer weights / (4*pi*r) with an
% explicit DiagonalValue for the r = 0 self entry.

properties (Constant)
    kernel = "single_layer_collocation_1_over_4pi_r"
    policy = "education_only_readable_hmatrix_not_production_solver"
end

properties
    targetPoints    % collocation points of the rows (nT x 3)
    sourcePoints    % collocation points of the columns (nS x 3)
    targetWeights   % row weights (areas for SurfaceMesh input)
    sourceWeights   % column quadrature weights
    targetTree      % row cluster tree
    sourceTree      % column cluster tree
    root            % block tree: dense / low_rank / split nodes
    leafSize        % cluster leaf size
    eta             % admissibility parameter
    rankTolerance   % relative SVD truncation tolerance
    diagonalValue   % kernel value used on the r = 0 diagonal
end

methods
    function H = HMatrix(target, source, options)
        arguments
            target
            source = []
            options.LeafSize (1,1) double {mustBeInteger, mustBePositive} = 8
            options.Eta (1,1) double {mustBePositive} = 2.0
            options.RankTolerance (1,1) double {mustBePositive} = 1e-6
            options.DiagonalValue (1,1) double = 0.0
        end

        [H.targetPoints, H.targetWeights] = bemCollocationPoints(target);
        if isempty(source)
            H.sourcePoints = H.targetPoints;
            H.sourceWeights = H.targetWeights;
        else
            [H.sourcePoints, H.sourceWeights] = bemCollocationPoints(source);
        end
        H.leafSize = options.LeafSize;
        H.eta = options.Eta;
        H.rankTolerance = options.RankTolerance;
        H.diagonalValue = options.DiagonalValue;

        H.targetTree = buildCluster(H.targetPoints, ...
            (1:size(H.targetPoints, 1)).', H.leafSize);
        H.sourceTree = buildCluster(H.sourcePoints, ...
            (1:size(H.sourcePoints, 1)).', H.leafSize);
        H.root = buildBlock(H.targetTree, H.sourceTree, ...
            H.targetPoints, H.sourcePoints, H.sourceWeights, ...
            H.eta, H.rankTolerance, H.diagonalValue);
    end

    function s = shape(H, dim)
        %SHAPE [nTarget, nSource] of the operator (not the object array size).
        s = [size(H.targetPoints, 1), size(H.sourcePoints, 1)];
        if nargin == 2
            s = s(dim);
        end
    end

    function y = matvec(H, x)
        %MATVEC Walk the block tree: dense, low-rank, and split cases.
        arguments
            H (1,1) HMatrix
            x (:,1) double
        end
        if numel(x) ~= H.shape(2)
            error("HMatrix:size", ...
                "Input vector length must match the H-matrix source size.");
        end
        y = zeros(H.shape(1), 1);
        y = applyBlock(H.root, x, y);
    end

    function y = mtimes(H, x)
        %MTIMES H * x is the block-tree matvec.
        y = H.matvec(x);
    end

    function stats = stats(H)
        %STATS Count blocks, ranks, and stored entries of the block tree.
        %
        % Use this to teach what compression did before worrying about
        % performance: dense near blocks stay dense, admissible far blocks
        % become low-rank, split blocks simply recurse.
        stats = struct();
        stats.blocks = 0;
        stats.denseBlocks = 0;
        stats.lowRankBlocks = 0;
        stats.splitBlocks = 0;
        stats.maxRank = 0;
        stats.storedEntries = 0;
        stats.compressionRatio = NaN;
        stats = visitBlock(H.root, stats);
        stats.compressionRatio = stats.storedEntries / prod(H.shape());
    end
end
end


function cluster = buildCluster(points, ids, leafSize)
%BUILDCLUSTER Recursive bisection cluster tree over point ids.

selected = points(ids, :);
cluster = struct();
cluster.ids = ids(:);
cluster.bboxMin = min(selected, [], 1);
cluster.bboxMax = max(selected, [], 1);
cluster.center = 0.5 * (cluster.bboxMin + cluster.bboxMax);
cluster.radius = max(sqrt(sum((selected - cluster.center).^2, 2)));
cluster.isLeaf = numel(ids) <= leafSize;
cluster.left = [];
cluster.right = [];

if cluster.isLeaf
    return
end

[~, axisId] = max(cluster.bboxMax - cluster.bboxMin);
[~, order] = sort(selected(:, axisId));
sortedIds = ids(order);
mid = floor(numel(sortedIds) / 2);
cluster.left = buildCluster(points, sortedIds(1:mid), leafSize);
cluster.right = buildCluster(points, sortedIds(mid + 1:end), leafSize);
end


function block = buildBlock(targetCluster, sourceCluster, targetPoints, sourcePoints, sourceWeights, eta, rankTolerance, diagonalValue)
%BUILDBLOCK Recursive block tree: low-rank far, dense near, split otherwise.

block = struct();
block.rows = targetCluster.ids(:);
block.cols = sourceCluster.ids(:);
block.children = {};

if isAdmissible(targetCluster, sourceCluster, eta)
    dense = laplaceBlock(targetPoints(block.rows, :), sourcePoints(block.cols, :), ...
        sourceWeights(block.cols), diagonalValue);
    [u, s, v] = svd(dense, "econ");
    singularValues = diag(s);
    if isempty(singularValues)
        rank = 0;
    else
        rank = find(singularValues > rankTolerance * singularValues(1), 1, "last");
    end
    if isempty(rank)
        rank = 0;
    end
    block.kind = "low_rank";
    block.rank = rank;
    block.U = u(:, 1:rank) * s(1:rank, 1:rank);
    block.V = v(:, 1:rank);
    return
end

if targetCluster.isLeaf && sourceCluster.isLeaf
    block.kind = "dense";
    block.rank = min(numel(block.rows), numel(block.cols));
    block.A = laplaceBlock(targetPoints(block.rows, :), sourcePoints(block.cols, :), ...
        sourceWeights(block.cols), diagonalValue);
    return
end

block.kind = "split";
block.rank = 0;
targetChildren = childClusters(targetCluster);
sourceChildren = childClusters(sourceCluster);
children = cell(numel(targetChildren), numel(sourceChildren));
for i = 1:numel(targetChildren)
    for j = 1:numel(sourceChildren)
        children{i, j} = buildBlock(targetChildren{i}, sourceChildren{j}, ...
            targetPoints, sourcePoints, sourceWeights, ...
            eta, rankTolerance, diagonalValue);
    end
end
block.children = children(:).';
end


function tf = isAdmissible(targetCluster, sourceCluster, eta)
%ISADMISSIBLE Far-field criterion: diameter <= eta * distance.

distance = norm(targetCluster.center - sourceCluster.center) ...
    - targetCluster.radius - sourceCluster.radius;
if distance <= 0
    tf = false;
    return
end
diameter = 2 * max(targetCluster.radius, sourceCluster.radius);
tf = diameter <= eta * distance;
end


function children = childClusters(cluster)
if cluster.isLeaf
    children = {cluster};
else
    children = {cluster.left, cluster.right};
end
end


function A = laplaceBlock(targetPoints, sourcePoints, sourceWeights, diagonalValue)
%LAPLACEBLOCK Dense collocation block weights / (4*pi*r).

nRows = size(targetPoints, 1);
nCols = size(sourcePoints, 1);
A = zeros(nRows, nCols);
for i = 1:nRows
    delta = sourcePoints - targetPoints(i, :);
    r = sqrt(sum(delta.^2, 2));
    values = sourceWeights ./ (4 * pi * r);
    values(r == 0) = diagonalValue;
    A(i, :) = values.';
end
end


function y = applyBlock(block, x, y)
%APPLYBLOCK Recursive matvec over dense / low_rank / split nodes.

switch block.kind
    case "dense"
        y(block.rows) = y(block.rows) + block.A * x(block.cols);
    case "low_rank"
        y(block.rows) = y(block.rows) + block.U * (block.V.' * x(block.cols));
    case "split"
        for k = 1:numel(block.children)
            y = applyBlock(block.children{k}, x, y);
        end
    otherwise
        error("HMatrix:block", "Unknown H-matrix block kind.");
end
end


function stats = visitBlock(block, stats)
%VISITBLOCK Accumulate per-kind block counts and stored entries.

stats.blocks = stats.blocks + 1;
switch block.kind
    case "dense"
        stats.denseBlocks = stats.denseBlocks + 1;
        stats.storedEntries = stats.storedEntries + numel(block.A);
    case "low_rank"
        stats.lowRankBlocks = stats.lowRankBlocks + 1;
        stats.maxRank = max(stats.maxRank, block.rank);
        stats.storedEntries = stats.storedEntries + numel(block.U) + numel(block.V);
    case "split"
        stats.splitBlocks = stats.splitBlocks + 1;
        for k = 1:numel(block.children)
            stats = visitBlock(block.children{k}, stats);
        end
    otherwise
        error("HMatrix:block", "Unknown H-matrix block kind.");
end
end
