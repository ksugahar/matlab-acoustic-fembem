function stats = educationalHMatrixStats(H)
%EDUCATIONALHMATRIXSTATS Count readable H-matrix block-tree statistics.
%
% Use this to teach what compression did before worrying about performance:
% dense near-field blocks stay dense, admissible far-field blocks become
% low-rank, and split blocks simply recurse.

arguments
    H (1,1) struct
end

stats = emptyStats();
stats = visitBlock(H.root, stats);
stats.compressionRatio = stats.storedEntries / prod(H.size);
end


function stats = emptyStats()
stats = struct();
stats.blocks = 0;
stats.denseBlocks = 0;
stats.lowRankBlocks = 0;
stats.splitBlocks = 0;
stats.maxRank = 0;
stats.storedEntries = 0;
stats.compressionRatio = NaN;
end


function stats = visitBlock(block, stats)
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
        error("educationalHMatrixStats:block", "Unknown H-matrix block kind.");
end
end
