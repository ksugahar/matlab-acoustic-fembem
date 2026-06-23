function H = educationalLaplaceHMatrix(target, source, options)
%EDUCATIONALLAPLACEHMATRIX Readable H-matrix for Laplace BEM education.
%
% This is intentionally small and inspectable. It demonstrates the core
% NGSolve.BEM/Gypsilab ideas: cluster trees, admissible far blocks, low-rank
% compression, dense near blocks, and matrix-free matvec.

arguments
    target
    source = []
    options.LeafSize (1,1) double {mustBeInteger, mustBePositive} = 8
    options.Eta (1,1) double {mustBePositive} = 2.0
    options.RankTolerance (1,1) double {mustBePositive} = 1e-6
    options.DiagonalValue (1,1) double = 0.0
end

[targetPoints, targetWeights] = extractBemPoints(target);
if isempty(source)
    sourcePoints = targetPoints;
    sourceWeights = targetWeights;
else
    [sourcePoints, sourceWeights] = extractBemPoints(source);
end

targetTree = buildCluster(targetPoints, (1:size(targetPoints, 1)).', options.LeafSize);
sourceTree = buildCluster(sourcePoints, (1:size(sourcePoints, 1)).', options.LeafSize);
root = buildBlock(targetTree, sourceTree, targetPoints, sourcePoints, sourceWeights, options);

H = struct();
H.kind = "educational_laplace_hmatrix";
H.kernel = "single_layer_collocation_1_over_4pi_r";
H.targetPoints = targetPoints;
H.sourcePoints = sourcePoints;
H.targetWeights = targetWeights;
H.sourceWeights = sourceWeights;
H.targetTree = targetTree;
H.sourceTree = sourceTree;
H.root = root;
H.options = options;
H.size = [size(targetPoints, 1), size(sourcePoints, 1)];
H.policy = "education_only_readable_hmatrix_not_production_solver";
end


function [points, weights] = extractBemPoints(input)
if isnumeric(input)
    points = input;
    weights = ones(size(points, 1), 1);
    return
end

if isstruct(input) && isfield(input, "gypsilab")
    vtx = input.gypsilab.vtx;
    tri = input.gypsilab.elt;
elseif isstruct(input) && isfield(input, "vtx") && isfield(input, "elt")
    vtx = input.vtx;
    tri = input.elt;
else
    error("educationalLaplaceHMatrix:input", ...
        "Input must be an Nx3 point array, a volFemBem model, or a struct with vtx/elt.");
end

points = (vtx(tri(:, 1), :) + vtx(tri(:, 2), :) + vtx(tri(:, 3), :)) / 3;
weights = triangleAreas(vtx, tri);
end


function areas = triangleAreas(vtx, tri)
a = vtx(tri(:, 1), :);
b = vtx(tri(:, 2), :);
c = vtx(tri(:, 3), :);
crossRows = cross(b - a, c - a, 2);
areas = 0.5 * sqrt(sum(crossRows.^2, 2));
end


function cluster = buildCluster(points, ids, leafSize)
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


function block = buildBlock(targetCluster, sourceCluster, targetPoints, sourcePoints, sourceWeights, options)
block = struct();
block.rows = targetCluster.ids(:);
block.cols = sourceCluster.ids(:);
block.children = {};

if isAdmissible(targetCluster, sourceCluster, options.Eta)
    dense = laplaceBlock(targetPoints(block.rows, :), sourcePoints(block.cols, :), ...
        sourceWeights(block.cols), options.DiagonalValue);
    [u, s, v] = svd(dense, "econ");
    singularValues = diag(s);
    if isempty(singularValues)
        rank = 0;
    else
        rank = find(singularValues > options.RankTolerance * singularValues(1), 1, "last");
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
        sourceWeights(block.cols), options.DiagonalValue);
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
            targetPoints, sourcePoints, sourceWeights, options);
    end
end
block.children = children(:).';
end


function tf = isAdmissible(targetCluster, sourceCluster, eta)
distance = norm(targetCluster.center - sourceCluster.center) - targetCluster.radius - sourceCluster.radius;
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
