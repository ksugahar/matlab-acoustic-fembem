function K = lowFrequencyStableHelmholtzKernel(targetPoints, sourcePoints, options)
%LOWFREQUENCYSTABLEHELMHOLTZKERNEL Stable dense Helmholtz BEM kernels.
%
% This educational helper mirrors the low-frequency-stable idea used by
% serious BEM implementations: keep the Laplace kernel explicit and evaluate
% the Helmholtz correction with expm1/Taylor series.
%
%   G_k(r) = exp(1i*k*r)/(4*pi*r)
%          = 1/(4*pi*r) + (exp(1i*k*r)-1)/(4*pi*r)
%
% For small k*r, the correction is evaluated as a power series so that the
% k -> 0 limit is smooth and readable.

arguments
    targetPoints (:,3) double
    sourcePoints (:,3) double
    options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
    options.SourceWeights double = []
    options.SourceNormals double = []
    options.DiagonalValue (1,1) double = 0.0
    options.DoubleLayerDiagonalValue (1,1) double = 0.0
    options.TaylorCutoff (1,1) double {mustBePositive} = 1e-4
    options.SeriesTerms (1,1) double {mustBeInteger, mustBePositive} = 8
end

nSource = size(sourcePoints, 1);
if isempty(options.SourceWeights)
    sourceWeights = ones(nSource, 1);
else
    sourceWeights = options.SourceWeights(:);
end
if numel(sourceWeights) ~= nSource
    error("lowFrequencyStableHelmholtzKernel:weights", ...
        "SourceWeights must have one entry per source point.");
end

hasNormals = ~isempty(options.SourceNormals);
if hasNormals && ~isequal(size(options.SourceNormals), size(sourcePoints))
    error("lowFrequencyStableHelmholtzKernel:normals", ...
        "SourceNormals must be an Nx3 array matching sourcePoints.");
end

nTarget = size(targetPoints, 1);
singleLaplace = zeros(nTarget, nSource);
singleCorrection = zeros(nTarget, nSource);
doubleLaplace = zeros(nTarget, nSource);
doubleCorrection = zeros(nTarget, nSource);

for i = 1:nTarget
    for j = 1:nSource
        delta = targetPoints(i, :) - sourcePoints(j, :);
        r = norm(delta);
        if r == 0
            singleLaplace(i, j) = options.DiagonalValue;
            if hasNormals
                doubleLaplace(i, j) = options.DoubleLayerDiagonalValue;
            end
            continue
        end

        weight = sourceWeights(j);
        z = 1i * options.Wavenumber * r;
        singleLaplace(i, j) = weight / (4 * pi * r);
        singleCorrection(i, j) = weight * stableExpm1OverR( ...
            z, r, options.TaylorCutoff, options.SeriesTerms) / (4 * pi);

        if hasNormals
            normalDot = dot(delta, options.SourceNormals(j, :));
            base = weight * normalDot / (4 * pi * r^3);
            factor = stableExpTimesOneMinusZMinusOne( ...
                z, options.TaylorCutoff, options.SeriesTerms);
            doubleLaplace(i, j) = base;
            doubleCorrection(i, j) = base * factor;
        end
    end
end

K = struct();
K.kind = "low_frequency_stable_helmholtz_kernel";
K.wavenumber = options.Wavenumber;
K.singleLayerLaplace = singleLaplace;
K.singleLayerCorrection = singleCorrection;
K.singleLayer = singleLaplace + singleCorrection;
K.policy = "low_frequency_stable_expm1_taylor_helmholtz_kernel";
if hasNormals
    K.doubleLayerSourceNormalLaplace = doubleLaplace;
    K.doubleLayerSourceNormalCorrection = doubleCorrection;
    K.doubleLayerSourceNormal = doubleLaplace + doubleCorrection;
end
end


function value = stableExpm1OverR(z, r, cutoff, terms)
if abs(z) < cutoff
    value = 0.0;
    for n = 1:terms
        value = value + z^n / factorial(n) / r;
    end
else
    value = expm1(z) / r;
end
end


function value = stableExpTimesOneMinusZMinusOne(z, cutoff, terms)
if abs(z) < cutoff
    value = 0.0;
    for n = 2:terms
        value = value + (1 - n) * z^n / factorial(n);
    end
else
    value = exp(z) * (1 - z) - 1;
end
end
