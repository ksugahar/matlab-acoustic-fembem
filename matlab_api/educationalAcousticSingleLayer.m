function op = educationalAcousticSingleLayer(target, source, options)
%EDUCATIONALACOUSTICSINGLELAYER Readable Helmholtz single-layer operator.
%
% This is the acoustic BEM teaching counterpart of the Laplace H-matrix demo.
% It keeps the kernel visible:
%
%   G_k(x,y) = exp(1i*k*|x-y|) / (4*pi*|x-y|)
%
% The implementation is dense and deliberately simple. Use it to understand
% NGSolve.BEM/Gypsilab acoustics before moving to compressed operators.

arguments
    target
    source = []
    options.Wavenumber (1,1) double {mustBeNonnegative} = 1.0
    options.DiagonalValue (1,1) double = 0.0
end

[targetPoints, targetWeights] = extractBemPoints(target);
if isempty(source)
    sourcePoints = targetPoints;
    sourceWeights = targetWeights;
else
    [sourcePoints, sourceWeights] = extractBemPoints(source);
end

kernelParts = lowFrequencyStableHelmholtzKernel(targetPoints, sourcePoints, ...
    "SourceWeights", sourceWeights, ...
    "Wavenumber", options.Wavenumber, ...
    "DiagonalValue", options.DiagonalValue);
A = kernelParts.singleLayer;

op = struct();
op.kind = "educational_acoustic_single_layer";
op.kernel = "exp(1i*k*r)/(4*pi*r)";
op.wavenumber = options.Wavenumber;
op.targetPoints = targetPoints;
op.sourcePoints = sourcePoints;
op.targetWeights = targetWeights;
op.sourceWeights = sourceWeights;
op.matrix = A;
op.kernelParts = kernelParts;
op.apply = @(x) A * x;
op.policy = "education_only_low_frequency_stable_dense_helmholtz_bem";
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
    error("educationalAcousticSingleLayer:input", ...
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
