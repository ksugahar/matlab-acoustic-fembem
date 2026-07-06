function report = structuredBoxVol(volFile, options)
%STRUCTUREDBOXVOL Structured rectangular-box tet mesh to .vol (no PDE Toolbox).
%
%   report = structuredBoxVol("box.vol", Size=[1.5 1.0 0.8], Cells=[6 4 4])
%
% Generates a first-order tetrahedral mesh of an axis-aligned box by splitting a
% Cells(1) x Cells(2) x Cells(3) grid of hexahedra into 6 conforming tetrahedra
% each (Kuhn decomposition around the local cell diagonal), then exports it with
% writePdeMeshVol.  Unlike writePdeBoxVol this needs NO MATLAB PDE Toolbox, so it
% runs anywhere the FEM/BEM lane runs -- handy for a rectangular scatterer whose
% whole boundary is the interface Gamma (the arbitrary-geometry counterpart to
% the analytic sphere fixtures).

arguments
    volFile (1,1) string
    options.Size (1,3) double {mustBePositive} = [1 1 1]
    options.Cells (1,3) double {mustBeInteger, mustBePositive} = [4 4 4]
    options.Centered (1,1) logical = true
    options.MaterialName (1,1) string = "domain"
    options.BoundaryName (1,1) string = "outer"
end

[nodes, tets] = boxTetMesh(options.Size, options.Cells);
if options.Centered
    nodes = nodes - mean([min(nodes); max(nodes)], 1);
end

report = writePdeMeshVol(struct("Nodes", nodes.', "Elements", tets.'), volFile, ...
    MaterialName=options.MaterialName, BoundaryName=options.BoundaryName);
report.generator = "structured_kuhn_6tet_box";
report.box_size = options.Size;
report.cells = options.Cells;
report.centered = options.Centered;
end


function [nodes, tets] = boxTetMesh(L, n)
%BOXTETMESH Grid of n hex cells over box L, each split into 6 conforming tets.
nx = n(1); ny = n(2); nz = n(3);
[X, Y, Z] = ndgrid(linspace(0,L(1),nx+1), linspace(0,L(2),ny+1), linspace(0,L(3),nz+1));
nodes = [X(:), Y(:), Z(:)];
gid = @(i,j,k) i + (nx+1)*((j-1) + (ny+1)*(k-1));   % 1-based, ndgrid(:) order
% 6-tet Kuhn split around the local (0,0,0)->(1,1,1) diagonal (corners v0..v7);
% the uniform global diagonal keeps shared inter-cell faces conforming.
splits = [1 2 3 7; 1 3 4 7; 1 4 8 7; 1 8 5 7; 1 5 6 7; 1 6 2 7];
tets = zeros(6*nx*ny*nz, 4);
e = 0;
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            c = [gid(i,j,k),   gid(i+1,j,k),   gid(i+1,j+1,k),   gid(i,j+1,k), ...
                 gid(i,j,k+1), gid(i+1,j,k+1), gid(i+1,j+1,k+1), gid(i,j+1,k+1)];
            for t = 1:6
                tet = c(splits(t, :));
                e1 = nodes(tet(2),:) - nodes(tet(1),:);
                e2 = nodes(tet(3),:) - nodes(tet(1),:);
                e3 = nodes(tet(4),:) - nodes(tet(1),:);
                if dot(cross(e1, e2), e3) < 0     % keep positive tet volume
                    tet([3 4]) = tet([4 3]);
                end
                e = e + 1;
                tets(e, :) = tet;
            end
        end
    end
end
end
