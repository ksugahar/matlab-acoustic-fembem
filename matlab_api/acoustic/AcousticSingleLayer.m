classdef AcousticSingleLayer
%ACOUSTICSINGLELAYER Readable dense Helmholtz single-layer operator.
%
%   op = AcousticSingleLayer(target, [], "Wavenumber", 10.0);
%   p = op.apply(q);         % or op * q
%   op.matrix                % the dense operator, kernel visible
%   op.kernelParts           % HelmholtzKernel with the Laplace/correction split
%
% This is the acoustic BEM teaching counterpart of the Laplace HMatrix demo:
%
%   G_k(x,y) = exp(1i*k*|x-y|) / (4*pi*|x-y|)
%
% assembled densely through the low-frequency-stable HelmholtzKernel so the
% k -> 0 limit stays connected to the Laplace BEM kernel. Understand the
% kernel and signs here, then compress with an H-matrix or move to
% NGSolve.BEM.

properties (Constant)
    kernel = "exp(1i*k*r)/(4*pi*r)"
    policy = "education_only_low_frequency_stable_dense_helmholtz_bem"
end

properties
    wavenumber      % k in the Helmholtz kernel
    targetPoints    % P1 nodal points of the rows (nT x 3)
    sourcePoints    % P1 nodal points of the columns (nS x 3)
    targetWeights   % row lumped P1 weights
    sourceWeights   % column quadrature weights
    matrix          % dense single-layer operator (nT x nS)
    kernelParts     % HelmholtzKernel: Laplace part + smooth correction
end

methods
    function op = AcousticSingleLayer(target, source, options)
        arguments
            target
            source = []
            options.Wavenumber (1,1) double {mustBeNonnegative} = 1.0
            options.DiagonalValue (1,1) double = 0.0
        end

        [op.targetPoints, op.targetWeights] = bemCollocationPoints(target);
        if isempty(source)
            op.sourcePoints = op.targetPoints;
            op.sourceWeights = op.targetWeights;
        else
            [op.sourcePoints, op.sourceWeights] = bemCollocationPoints(source);
        end

        op.wavenumber = options.Wavenumber;
        op.kernelParts = HelmholtzKernel(op.targetPoints, op.sourcePoints, ...
            "SourceWeights", op.sourceWeights, ...
            "Wavenumber", options.Wavenumber, ...
            "DiagonalValue", options.DiagonalValue);
        op.matrix = op.kernelParts.singleLayer;
    end

    function s = shape(op, dim)
        %SHAPE [nTarget, nSource] of the operator (not the object array size).
        s = size(op.matrix);
        if nargin == 2
            s = s(dim);
        end
    end

    function p = apply(op, q)
        %APPLY Dense single-layer potential of a source density.
        p = op.matrix * q;
    end

    function p = mtimes(op, q)
        %MTIMES op * q is the same dense apply.
        p = op.matrix * q;
    end
end
end
