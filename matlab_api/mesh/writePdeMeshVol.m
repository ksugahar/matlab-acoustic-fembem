function report = writePdeMeshVol(pdeMesh, volFile, options)
%WRITEPDEMESHVOL Export a MATLAB PDE Toolbox tetrahedral mesh to .vol.
%
%   model = createpde();
%   model.Geometry = multicuboid(1, 1, 1);
%   generateMesh(model, "GeometricOrder", "linear");
%   report = writePdeMeshVol(model.Mesh, "box.vol");
%
% The FEM/BEM teaching lane accepts first-order Netgen .vol files with
% triangle boundary faces and tetrahedral volume cells.  This exporter keeps
% that contract strict: PDE Toolbox meshes must be linear tet meshes.

arguments
    pdeMesh
    volFile (1,1) string
    options.MaterialId (1,1) double {mustBeInteger, mustBePositive} = 1
    options.MaterialName (1,1) string = "domain"
    options.BoundaryId (1,1) double {mustBeInteger, mustBePositive} = 1
    options.BoundaryName (1,1) string = "outer"
end

[points, tets] = pdeNodesAndTets(pdeMesh);
tri = boundaryTriangles(points, tets);

% Delegate the .vol formatting to the pure-MATLAB array writer so there is a
% single canonical writer; this adapter only turns a PDE mesh into arrays.
bnames = strings(1, options.BoundaryId);
bnames(options.BoundaryId) = options.BoundaryName;
mnames = strings(1, options.MaterialId);
mnames(options.MaterialId) = options.MaterialName;
wr = writeVol(volFile, points, tri, ...
    Tets=tets, ...
    TriBoundaryId=repmat(options.BoundaryId, size(tri, 1), 1), ...
    TetMaterialId=repmat(options.MaterialId, size(tets, 1), 1), ...
    BoundaryNames=bnames, ...
    MaterialNames=mnames);
report = struct();
report.tool = "write_pde_mesh_vol";
report.status = "ok";
report.policy = "pde_toolbox_linear_tet_to_netgen_vol_tri_tet";
report.output_file = volFile;
report.input_family = "matlab_pde_toolbox_mesh";
report.points = size(points, 1);
report.triangles = size(tri, 1);
report.tets = size(tets, 1);
report.material_id = options.MaterialId;
report.material_name = options.MaterialName;
report.boundary_id = options.BoundaryId;
report.boundary_name = options.BoundaryName;
report.roundtrip_summary = wr.roundtrip_summary;
report.boundary_orientation = wr.boundary_orientation;
end


function [points, tets] = pdeNodesAndTets(pdeMesh)
if isobject(pdeMesh) || isstruct(pdeMesh)
    nodes = pdeMesh.Nodes;
    elements = pdeMesh.Elements;
else
    error("writePdeMeshVol:input", ...
        "pdeMesh must be a PDE Toolbox mesh object or a struct with Nodes and Elements.");
end

if size(nodes, 1) == 3
    points = double(nodes).';
elseif size(nodes, 2) == 3
    points = double(nodes);
else
    error("writePdeMeshVol:nodes", ...
        "PDE mesh Nodes must be 3-by-N or N-by-3 coordinates.");
end

if size(elements, 1) == 4
    tets = double(elements).';
elseif size(elements, 2) == 4
    tets = double(elements);
else
    error("writePdeMeshVol:elementOrder", ...
        "Only linear tetrahedral PDE meshes are accepted. Use generateMesh(..., ""GeometricOrder"", ""linear"").");
end

if any(tets(:) < 1) || any(tets(:) > size(points, 1))
    error("writePdeMeshVol:nodes", "Element connectivity references nodes outside 1..N.");
end
end


function tri = boundaryTriangles(points, tets)
localFaces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
nTet = size(tets, 1);
faces = zeros(4 * nTet, 3);
adjacentTet = zeros(4 * nTet, 1);
row = 0;
for k = 1:nTet
    for f = 1:4
        row = row + 1;
        faces(row, :) = tets(k, localFaces(f, :));
        adjacentTet(row) = k;
    end
end

[~, ~, faceGroup] = unique(sort(faces, 2), "rows");
faceCounts = accumarray(faceGroup, 1);
isBoundary = faceCounts(faceGroup) == 1;
tri = faces(isBoundary, :);
adjacentTet = adjacentTet(isBoundary);

for k = 1:size(tri, 1)
    face = tri(k, :);
    tet = tets(adjacentTet(k), :);
    areaVec = 0.5 * cross( ...
        points(face(2), :) - points(face(1), :), ...
        points(face(3), :) - points(face(1), :));
    faceCentroid = mean(points(face, :), 1);
    tetCentroid = mean(points(tet, :), 1);
    if dot(areaVec, faceCentroid - tetCentroid) < 0
        tri(k, [2 3]) = tri(k, [3 2]);
    end
end
end
