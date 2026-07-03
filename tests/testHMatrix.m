function tests = testHMatrix
%TESTHMATRIX Tests for readable H-matrix scaffolds.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testSeparatedClustersMatvecMatchesDenseKernel(testCase)
target = [ ...
    0.0 0.0 0.0
    0.2 0.0 0.0
    0.4 0.0 0.0
    0.6 0.0 0.0];
source = target + [3.0 0.5 0.0];

H = HMatrix(target, source, ...
    "LeafSize", 2, "Eta", 1.0, "RankTolerance", 1e-14);
x = (1:4).';

verifyEqual(testCase, H.matvec(x), directLaplace(target, source) * x, ...
    "AbsTol", 1e-12);
verifyError(testCase, @() H.matvec(ones(5, 1)), "HMatrix:size");

stats = H.stats();
verifyGreaterThan(testCase, stats.lowRankBlocks, 0);
verifyEqual(testCase, H.policy, "education_only_readable_hmatrix_not_production_solver");
end


function testSurfaceP1NodesBecomeWeightedBemPoints(testCase)
path = writeFixture(testCase, tetVolText());
mesh = VolMesh(path);
surface = mesh.boundary();
expectedWeights = [0.5; repmat((1 + sqrt(3) / 2) / 3, 3, 1)];

H = HMatrix(surface, [], ...
    "LeafSize", 2, "Eta", 1.5, "RankTolerance", 1e-12);
y = H.matvec(ones(4, 1));
stats = H.stats();

verifyEqual(testCase, H.shape(), [4 4]);
verifyEqual(testCase, H.shape(2), 4);
verifyEqual(testCase, H.leafSize, 2);
verifyEqual(testCase, H.eta, 1.5);
verifyEqual(testCase, H.rankTolerance, 1e-12);
verifySize(testCase, H.targetPoints, [4 3]);
verifySize(testCase, H.sourceWeights, [4 1]);
verifyEqual(testCase, H.sourceWeights, expectedWeights, "AbsTol", 1e-14);
verifyEqual(testCase, H.kernel, "single_layer_p1_nodal_lumped_1_over_4pi_r");
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


function path = writeFixture(testCase, text)
path = string(fullfile(tempdir, "hmatrixTest_" + char(java.util.UUID.randomUUID()) + ".vol"));
fid = fopen(path, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", text);
clear cleanup
testCase.addTeardown(@() delete(path));
end


function text = tetVolText()
text = join([
    "mesh3d"
    "dimension"
    "3"
    "geomtype"
    "0"
    "facedescriptors"
    "1"
    "1 1 0 1 1"
    "surfaceelements"
    "4"
    "1 1 1 0 3 1 2 3"
    "1 1 1 0 3 1 4 2"
    "1 1 1 0 3 2 4 3"
    "1 1 1 0 3 3 4 1"
    "volumeelements"
    "1"
    "1 4 1 2 3 4"
    "points"
    "4"
    "0 0 0"
    "1 0 0"
    "0 1 0"
    "0 0 1"
    "pointelements"
    "0"
    "materials"
    "1"
    "1 air"
    "bcnames"
    "1"
    "1 outer"
    "endmesh"
    ], newline);
end
