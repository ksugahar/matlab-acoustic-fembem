function study = hmatrixScalingStudy(pointCounts, options)
%HMATRIXSCALINGSTUDY Compare readable ACA+ H-matrices with dense kernels.
arguments
    pointCounts (1,:) double {mustBeInteger, mustBePositive} = [60 120 240]
    options.LeafSize (1,1) double {mustBeInteger, mustBePositive} = 32
    options.Eta (1,1) double {mustBePositive} = 2.0
    options.RankTolerance (1,1) double {mustBePositive} = 1e-8
end
if numel(pointCounts) < 3 || any(diff(pointCounts) <= 0)
    error("hmatrixScalingStudy:counts", ...
        "Use at least three strictly increasing point counts.");
end

rows = repmat(struct(), numel(pointCounts), 1);
for k = 1:numel(pointCounts)
    n = pointCounts(k);
    t = linspace(0, 1, n).';
    target = [t, zeros(n, 1), zeros(n, 1)];
    source = target + [5, 0.25, 0.1];
    x = sin((1:n).' * sqrt(2)) + 0.25*cos((1:n).' * sqrt(3));

    started = tic;
    H = HMatrix(target, source, ...
        "LeafSize", options.LeafSize, ...
        "Eta", options.Eta, ...
        "RankTolerance", options.RankTolerance);
    rows(k).buildSeconds = toc(started);

    started = tic;
    yH = H.matvec(x);
    rows(k).matvecSeconds = toc(started);

    started = tic;
    A = directLaplace(target, source);
    yDense = A*x;
    rows(k).denseReferenceSeconds = toc(started);

    stats = H.stats();
    rows(k).pointCount = n;
    rows(k).maxRank = stats.maxRank;
    rows(k).lowRankBlocks = stats.lowRankBlocks;
    rows(k).denseBlocks = stats.denseBlocks;
    rows(k).storedEntries = stats.storedEntries;
    rows(k).denseEntries = n*n;
    rows(k).compressionRatio = stats.compressionRatio;
    rows(k).matvecRelativeError = norm(yH-yDense)/norm(yDense);
end

study = struct( ...
    "schema", "matlab-acoustic-fembem.hmatrix-scaling-study.v1", ...
    "kernel", HMatrix.kernel, ...
    "policy", HMatrix.policy, ...
    "leafSize", options.LeafSize, ...
    "eta", options.Eta, ...
    "rankTolerance", options.RankTolerance, ...
    "rows", rows);
end


function A = directLaplace(target, source)
nRows = size(target, 1);
nCols = size(source, 1);
A = zeros(nRows, nCols);
for i = 1:nRows
    delta = source - target(i, :);
    r = sqrt(sum(delta.^2, 2));
    A(i, :) = (1 ./ (4*pi*r)).';
end
end
