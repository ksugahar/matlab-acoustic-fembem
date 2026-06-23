function y = educationalHMatrixMatvec(H, x)
%EDUCATIONALHMATRIXMATVEC Apply the readable educational H-matrix.
%
% The implementation favors clarity over speed. It walks the block tree and
% makes the dense/low-rank/split cases explicit, mirroring the ideas behind
% Gypsilab hmx and NGSolve.BEM compression.

arguments
    H (1,1) struct
    x (:,1) double
end

if numel(x) ~= H.size(2)
    error("educationalHMatrixMatvec:size", ...
        "Input vector length must match the H-matrix source size.");
end

y = zeros(H.size(1), 1);
y = applyBlock(H.root, x, y);
end


function y = applyBlock(block, x, y)
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
        error("educationalHMatrixMatvec:block", "Unknown H-matrix block kind.");
end
end
