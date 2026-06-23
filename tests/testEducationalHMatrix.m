function tests = testEducationalHMatrix
%TESTEDUCATIONALHMATRIX Tests for readable educational H-matrix scaffolds.

tests = functiontests(localfunctions);
end


function testSeparatedClustersMatvecMatchesDenseKernel(testCase)
target = [ ...
    0.0 0.0 0.0
    0.2 0.0 0.0
    0.4 0.0 0.0
    0.6 0.0 0.0];
source = target + [3.0 0.5 0.0];

H = educationalLaplaceHMatrix(target, source, ...
    "LeafSize", 2, "Eta", 1.0, "RankTolerance", 1e-14);
x = (1:4).';

verifyEqual(testCase, educationalHMatrixMatvec(H, x), directLaplace(target, source) * x, ...
    "AbsTol", 1e-12);

stats = educationalHMatrixStats(H);
verifyGreaterThan(testCase, stats.lowRankBlocks, 0);
verifyEqual(testCase, H.policy, "education_only_readable_hmatrix_not_production_solver");
end


function testSurfaceTrianglesBecomeWeightedBemPoints(testCase)
surface = struct();
surface.vtx = [ ...
    0 0 0
    1 0 0
    0 1 0
    0 0 1];
surface.elt = [ ...
    1 2 3
    1 4 2
    2 4 3
    3 4 1];

H = educationalLaplaceHMatrix(surface, [], ...
    "LeafSize", 2, "Eta", 1.5, "RankTolerance", 1e-12);
y = educationalHMatrixMatvec(H, ones(4, 1));
stats = educationalHMatrixStats(H);

verifySize(testCase, H.targetPoints, [4 3]);
verifySize(testCase, H.sourceWeights, [4 1]);
verifyEqual(testCase, H.sourceWeights(1:2), [0.5; 0.5], "AbsTol", 1e-14);
verifySize(testCase, y, [4 1]);
verifyGreaterThan(testCase, stats.blocks, 0);
verifyLessThanOrEqual(testCase, stats.compressionRatio, 1.0);
end


function A = directLaplace(target, source)
A = zeros(size(target, 1), size(source, 1));
for i = 1:size(target, 1)
    delta = source - target(i, :);
    r = sqrt(sum(delta.^2, 2));
    A(i, :) = (1 ./ (4 * pi * r)).';
end
end
