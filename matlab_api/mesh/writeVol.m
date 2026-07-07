function report = writeVol(volFile, points, tri, options)
%WRITEVOL Write a first-order Netgen .vol from raw arrays (pure MATLAB).
%
%   writeVol("m.vol", points, tri)                 % surface-only (BEM-ready)
%   writeVol("m.vol", points, tri, Tets=tets)      % + tetrahedra (FEM volume)
%   writeVol("m.vol", points, tri, ...
%       TriBoundaryId=col, BoundaryNames=["inner" "outer"], ...
%       TetMaterialId=mat, MaterialNames=["iron" "air"])
%
% This is the WRITE counterpart of readVolTriTet / VolMesh: MATLAB owns .vol
% read AND write with no PDE Toolbox and no netgen.  It emits the same
% first-order "mesh3d" .vol the readers accept (triangle boundary faces +
% optional tetrahedra), so a mesh built by ANY means -- a pure-MATLAB
% tessellation, an edited mesh, an analytic shape -- can be written to .vol.
%
%   points        nNodes x 3 node coordinates
%   tri           nTri x 3 one-based boundary-triangle node ids
%   Tets          nTet x 4 one-based tetra node ids (default none -> surface-only)
%   TriBoundaryId nTri x 1 Netgen boundary number per triangle (default all 1)
%   TetMaterialId nTet x 1 Netgen material number per tet     (default all 1)
%   BoundaryNames string array; BoundaryNames(id) names boundary id (default "outer")
%   MaterialNames string array; MaterialNames(id) names material id (default "domain")
%
% Surface-only .vol (no Tets) is valid: SurfaceMesh / the BEM spaces read the
% boundary triangles directly and never need the interior tets.  The written
% file is round-tripped through readVolTriTet before returning, so a malformed
% write fails loudly here rather than at the next read.

arguments
    volFile (1,1) string
    points (:,3) double {mustBeFinite}
    tri (:,3) double {mustBeInteger, mustBePositive}
    options.Tets (:,4) double {mustBeInteger, mustBePositive} = zeros(0, 4)
    options.TriBoundaryId (:,1) double {mustBeInteger, mustBePositive} = []
    options.TetMaterialId (:,1) double {mustBeInteger, mustBePositive} = []
    options.BoundaryNames (1,:) string = "outer"
    options.MaterialNames (1,:) string = "domain"
end

tets = options.Tets;
triCol = options.TriBoundaryId;
if isempty(triCol), triCol = ones(size(tri, 1), 1); end
tetMat = options.TetMaterialId;
if isempty(tetMat), tetMat = ones(size(tets, 1), 1); end
if numel(triCol) ~= size(tri, 1)
    error("writeVol:triId", "TriBoundaryId needs one id per triangle (%d).", size(tri, 1));
end
if numel(tetMat) ~= size(tets, 1)
    error("writeVol:tetId", "TetMaterialId needs one id per tetrahedron (%d).", size(tets, 1));
end
nNodes = size(points, 1);
if any(tri(:) > nNodes) || (~isempty(tets) && any(tets(:) > nNodes))
    error("writeVol:nodes", "Connectivity references nodes outside 1..%d.", nNodes);
end

outDir = fileparts(volFile);
if outDir ~= "" && ~isfolder(outDir)
    mkdir(outDir);
end
fid = fopen(volFile, "w");
if fid < 0
    error("writeVol:file", "Cannot open .vol output: %s", volFile);
end
cleanup = onCleanup(@() fclose(fid));

bcIds = unique(triCol(:)).';
matIds = unique([tetMat(:); 1]).';                 % always name at least material 1

fprintf(fid, "mesh3d\n");
fprintf(fid, "dimension\n3\n");
fprintf(fid, "geomtype\n0\n");
fprintf(fid, "facedescriptors\n%d\n", numel(bcIds));
for b = bcIds
    fprintf(fid, "%d %d 0 1 1\n", b, b);
end
fprintf(fid, "surfaceelements\n%d\n", size(tri, 1));
for k = 1:size(tri, 1)
    fprintf(fid, "1 %d 1 0 3 %d %d %d\n", triCol(k), tri(k, :));
end
fprintf(fid, "volumeelements\n%d\n", size(tets, 1));
for k = 1:size(tets, 1)
    fprintf(fid, "%d 4 %d %d %d %d\n", tetMat(k), tets(k, :));
end
fprintf(fid, "points\n%d\n", nNodes);
for k = 1:nNodes
    fprintf(fid, "%.17g %.17g %.17g\n", points(k, :));
end
fprintf(fid, "pointelements\n0\n");
fprintf(fid, "materials\n%d\n", numel(matIds));
for m = matIds
    fprintf(fid, "%d %s\n", m, char(nameFor(options.MaterialNames, m, "domain")));
end
fprintf(fid, "bcnames\n%d\n", numel(bcIds));
for b = bcIds
    fprintf(fid, "%d %s\n", b, char(nameFor(options.BoundaryNames, b, "outer")));
end
fprintf(fid, "endmesh\n");
clear cleanup

roundtrip = readVolTriTet(volFile);                % fail loudly on a malformed write
report = struct();
report.tool = "write_vol";
report.status = "ok";
report.policy = "pure_matlab_arrays_to_netgen_vol_tri_tet";
report.output_file = volFile;
report.points = nNodes;
report.triangles = size(tri, 1);
report.tets = size(tets, 1);
report.boundaries = numel(bcIds);
report.materials = numel(matIds);
report.surface_only = isempty(tets);
report.roundtrip_summary = roundtrip.summary;
report.boundary_orientation = roundtrip.boundaryOrientation.boundaryOrientation;
end


function nm = nameFor(names, id, fallback)
%NAMEFOR Name for a Netgen id from a 1-based string array, else a fallback.
if id >= 1 && id <= numel(names) && strlength(names(id)) > 0
    nm = names(id);
else
    nm = fallback + "_" + id;
end
end
