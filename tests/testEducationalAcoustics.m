function tests = testEducationalAcoustics
%TESTEDUCATIONALACOUSTICS Tests for readable acoustic BEM scaffolds.

tests = functiontests(localfunctions);
end


function testHelmholtzKernelMatchesManualTwoPointValue(testCase)
target = [0 0 0];
source = [2 0 0];
k = 3.0;

op = educationalAcousticSingleLayer(target, source, "Wavenumber", k);
expected = exp(1i * k * 2) / (4 * pi * 2);

verifyEqual(testCase, op.matrix, expected, "AbsTol", 1e-14);
verifyEqual(testCase, op.apply(2), 2 * expected, "AbsTol", 1e-14);
verifyEqual(testCase, op.policy, "education_only_low_frequency_stable_dense_helmholtz_bem");
end


function testLowFrequencyCorrectionMatchesTaylorLimit(testCase)
target = [0 0 0];
source = [0.75 0 0];
k = 1e-9;

K = lowFrequencyStableHelmholtzKernel(target, source, "Wavenumber", k);
expectedCorrection = 1i * k / (4 * pi) - k^2 * 0.75 / (8 * pi);

verifyEqual(testCase, K.singleLayerLaplace, 1 / (4 * pi * 0.75), "AbsTol", 1e-14);
verifyEqual(testCase, K.singleLayerCorrection, expectedCorrection, "RelTol", 1e-12);
verifyEqual(testCase, K.policy, "low_frequency_stable_expm1_taylor_helmholtz_kernel");
end


function testLowFrequencyDoubleLayerMatchesDirectFormula(testCase)
target = [0 0 0];
source = [2 0 0];
normal = [-1 0 0];
k = 1e-7;

K = lowFrequencyStableHelmholtzKernel(target, source, ...
    "Wavenumber", k, "SourceNormals", normal);

delta = target - source;
r = norm(delta);
expected = dot(delta, normal) * exp(1i * k * r) * (1 - 1i * k * r) / (4 * pi * r^3);

verifyEqual(testCase, K.doubleLayerSourceNormal, expected, "RelTol", 1e-12);
verifyEqual(testCase, K.doubleLayerSourceNormalCorrection, 0, "AbsTol", 1e-12);
end


function testZeroWavenumberReducesToLaplaceKernel(testCase)
target = [0 0 0; 0 1 0];
source = [2 0 0; 2 1 0];

op = educationalAcousticSingleLayer(target, source, "Wavenumber", 0.0);
expected = directLaplace(target, source);

verifyEqual(testCase, op.matrix, expected, "AbsTol", 1e-14);
end


function testSurfaceTrianglesUseAreaWeights(testCase)
surface = struct();
surface.vtx = [ ...
    0 0 0
    1 0 0
    0 1 0];
surface.elt = [1 2 3];

op = educationalAcousticSingleLayer(surface, surface, ...
    "Wavenumber", 0.0, "DiagonalValue", 0.25);

verifyEqual(testCase, op.sourceWeights, 0.5, "AbsTol", 1e-14);
verifyEqual(testCase, op.matrix, 0.25, "AbsTol", 1e-14);
end


function A = directLaplace(target, source)
A = zeros(size(target, 1), size(source, 1));
for i = 1:size(target, 1)
    delta = source - target(i, :);
    r = sqrt(sum(delta.^2, 2));
    A(i, :) = (1 ./ (4 * pi * r)).';
end
end
