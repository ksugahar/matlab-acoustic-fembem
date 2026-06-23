function model = volFemBemModel(volFile, gypsilabRoot)
%VOLFEMBEMMODEL Build a readable Lukas FEM / Gypsilab BEM mesh model.
%
% This function performs mesh intake only. It keeps the FEM/BEM trace explicit
% so scalar and later H(curl)/RWG coupling code can be tested independently.

arguments
    volFile (1,1) string
    gypsilabRoot (1,1) string = ""
end

if gypsilabRoot ~= ""
    addpath(fullfile(gypsilabRoot, "openMsh"));
    addpath(fullfile(gypsilabRoot, "openDom"));
    addpath(fullfile(gypsilabRoot, "openFem"));
    addpath(fullfile(gypsilabRoot, "openOpr"));
    addpath(fullfile(gypsilabRoot, "openHmx"));
end

mesh = readVolTriTet(volFile);
[boundaryVtx, boundaryTri, traceNodeIds] = compactBoundaryMesh(mesh);

model = struct();
model.mesh = mesh;
model.lukas.geo = makeLukasGeo(mesh);
model.gypsilab.vtx = boundaryVtx;
model.gypsilab.elt = boundaryTri;
model.gypsilab.col = mesh.triCol;
model.gypsilab.boundaryMesh = [];
model.gypsilab.globalNodeIds = traceNodeIds;
model.trace.nodeIds = traceNodeIds;
model.trace.surfaceTrianglesGlobal = mesh.tri;
model.trace.surfaceTrianglesLocal = boundaryTri;
model.trace.boundaryNumbers = mesh.triCol;
model.trace.policy = mesh.policy;
model.status = "mesh_ready";

if exist("msh", "file") == 2
    model.gypsilab.boundaryMesh = msh(boundaryVtx, boundaryTri, mesh.triCol);
    model.status = "mesh_ready_gypsilab_msh_constructed";
end
end


function [boundaryVtx, boundaryTri, traceNodeIds] = compactBoundaryMesh(mesh)
%COMPACTBOUNDARYMESH Keep Gypsilab P1 dofs compact and map back to Lukas nodes.

traceNodeIds = unique(mesh.tri(:));
globalToLocal = zeros(max(traceNodeIds), 1);
globalToLocal(traceNodeIds) = 1:numel(traceNodeIds);
boundaryTri = globalToLocal(mesh.tri);
boundaryVtx = mesh.vtx(traceNodeIds, :);
end


function geo = makeLukasGeo(mesh)
%MAKELUKASGEO Build the subset of Lukas FEM geo fields needed by prototypes.

geo = struct();
geo.nodes = mesh.vtx;
geo.conn_matrix = mesh.tet;
geo.M = size(mesh.tet, 1);
geo.N = size(mesh.vtx, 1);
geo.x = mesh.vtx(:, 1);
geo.y = mesh.vtx(:, 2);
geo.z = mesh.vtx(:, 3);
geo.regions = makeRegions(mesh);
geo.policy.order = 1;
geo.policy.family = "P1 tetrahedra";
end


function regions = makeRegions(mesh)
%MAKEREGIONS Group tetrahedra by material number.

matIds = unique(mesh.tetMat(:)).';
regions = struct("matnr", {}, "name", {}, "Elements", {}, "Nodes", {});
for k = 1:numel(matIds)
    matnr = matIds(k);
    elems = find(mesh.tetMat == matnr);
    nodes = unique(mesh.tet(elems, :));
    regions(k).matnr = matnr;
    if isKey(mesh.materials, matnr)
        regions(k).name = string(mesh.materials(matnr));
    else
        regions(k).name = "material_" + matnr;
    end
    regions(k).Elements = elems(:);
    regions(k).Nodes = nodes(:);
end
end
