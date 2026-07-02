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
boundaryNames = boundaryNamesFor(mesh.triCol, mesh.boundaries);
boundaryRowIdentity = makeBoundaryRowIdentity( ...
    mesh.tri, ...
    mesh.triCol, ...
    boundaryNames, ...
    mesh.boundaryOrientation.adjacentTetIndices);
sourcePath = string(volFile);
sourceFileId = localFileSha256Id(sourcePath);
[~, sourceStem, sourceExt] = fileparts(sourcePath);
meshId = "netgen_vol:" + string(sourceStem) + string(sourceExt);
surfaceMeshId = meshId + ":boundary_tri_p1";
traceArtifactId = meshId + ":h1_to_scalar_bem_trace_p1";
traceOperatorArtifactId = meshId + ":h1_to_scalar_bem_trace_operator_p1";
traceOperatorPolicy = "one_hot_boundary_node_injection_from_vol_node_ids";
traceOutputArtifactId = meshId + ":h1_to_scalar_bem_trace_output_p1";
traceOutputPath = "memory://" + traceOutputArtifactId;
traceObservableId = meshId + ":h1_to_scalar_bem_boundary_trace_observable_p1";
traceObservableFamily = "fem_bem_boundary_trace";
traceBasisSchemaId = "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1";
assemblyRuleId = "first_order_tet_h1_trace_tri_p1_bem_teaching_v1";
quadratureRuleId = "tri_p1_exact_mass_regular_kernel_teaching_v1";

model = struct();
model.identity = struct( ...
    "meshId", meshId, ...
    "surfaceMeshId", surfaceMeshId, ...
    "traceArtifactId", traceArtifactId, ...
    "traceOperatorArtifactId", traceOperatorArtifactId, ...
    "traceOperatorPolicy", traceOperatorPolicy, ...
    "traceOutputArtifactId", traceOutputArtifactId, ...
    "traceOutputPath", traceOutputPath, ...
    "traceObservableId", traceObservableId, ...
    "traceObservableFamily", traceObservableFamily, ...
    "traceBasisSchemaId", traceBasisSchemaId, ...
    "assemblyRuleId", assemblyRuleId, ...
    "quadratureRuleId", quadratureRuleId, ...
    "sourcePath", sourcePath, ...
    "sourceFileId", sourceFileId, ...
    "sourceFormat", ".vol");
model.mesh = mesh;
model.lukas.geo = makeLukasGeo(mesh);
model.gypsilab.vtx = boundaryVtx;
model.gypsilab.elt = boundaryTri;
model.gypsilab.col = mesh.triCol;
model.gypsilab.boundaryNumbers = mesh.triCol;
model.gypsilab.boundaryNames = boundaryNames;
model.gypsilab.boundaryRowIdentity = boundaryRowIdentity;
model.gypsilab.boundaryMesh = [];
model.gypsilab.globalNodeIds = traceNodeIds;
model.trace.nodeIds = traceNodeIds;
model.trace.meshId = meshId;
model.trace.surfaceMeshId = surfaceMeshId;
model.trace.traceArtifactId = traceArtifactId;
model.trace.traceOperatorArtifactId = traceOperatorArtifactId;
model.trace.traceOperatorPolicy = traceOperatorPolicy;
model.trace.traceOutputArtifactId = traceOutputArtifactId;
model.trace.traceOutputDigest = "";
model.trace.traceOutputPath = traceOutputPath;
model.trace.traceObservableId = traceObservableId;
model.trace.traceObservableFamily = traceObservableFamily;
model.trace.traceBasisSchemaId = traceBasisSchemaId;
model.trace.assemblyRuleId = assemblyRuleId;
model.trace.quadratureRuleId = quadratureRuleId;
model.trace.sourcePath = sourcePath;
model.trace.sourceFileId = sourceFileId;
model.trace.sourceFormat = ".vol";
model.trace.surfaceTrianglesGlobal = mesh.tri;
model.trace.surfaceTrianglesLocal = boundaryTri;
model.trace.traceRowIdentity = makeTraceRowIdentity(traceNodeIds);
model.trace.boundaryNumbers = mesh.triCol;
model.trace.boundaryNames = boundaryNames;
model.trace.boundaryRowIdentity = boundaryRowIdentity;
model.trace.boundaryOrientation = mesh.boundaryOrientation.boundaryOrientation;
model.trace.triangleOrientationSignsToOutward = mesh.boundaryOrientation.triangleOrientationSignsToOutward;
model.trace.adjacentTetIndices = mesh.boundaryOrientation.adjacentTetIndices;
model.trace.orientationRows = mesh.boundaryOrientation.rows;
model.trace.policy = mesh.policy;
model.status = "mesh_ready";

if exist("msh", "file") == 2
    model.gypsilab.boundaryMesh = msh(boundaryVtx, boundaryTri, mesh.triCol);
    model.status = "mesh_ready_gypsilab_msh_constructed";
end


function names = boundaryNamesFor(boundaryNumbers, boundaryMap)
%BOUNDARYNAMESFOR Expand Netgen boundary ids to per-triangle labels.

names = strings(numel(boundaryNumbers), 1);
for k = 1:numel(boundaryNumbers)
    id = boundaryNumbers(k);
    if isKey(boundaryMap, id)
        names(k) = string(boundaryMap(id));
    else
        names(k) = "boundary_" + id;
    end
end
end


function id = localFileSha256Id(path)
%LOCALFILESHA256ID Stable source identity for the .vol-to-trace handoff.

fid = fopen(path, "rb");
if fid < 0
    error("volFemBemModel:file", "Cannot open .vol source: %s", path);
end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, "*uint8");
md = javaMethod("getInstance", "java.security.MessageDigest", "SHA-256");
md.update(typecast(bytes(:), "int8"));
hash = typecast(md.digest(), "uint8");
hex = lower(reshape(dec2hex(hash, 2).', 1, []));
id = "sha256:" + string(hex);
clear cleanup
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


function identity = makeTraceRowIdentity(traceNodeIds)
%MAKETRACEROWIDENTITY Bind each compact boundary row to its FEM/BEM node ids.

nRows = numel(traceNodeIds);
identity = struct( ...
    "trace_row_index", num2cell((1:nRows).'), ...
    "fem_node_id", num2cell(traceNodeIds(:)), ...
    "bem_node_id", num2cell(traceNodeIds(:)), ...
    "surface_node_index", num2cell((1:nRows).'));
end


function identity = makeBoundaryRowIdentity(surfaceTriangles, boundaryNumbers, boundaryNames, adjacentTetIndices)
%MAKEBOUNDARYROWIDENTITY Bind each boundary triangle row to its BC labels.

nRows = size(surfaceTriangles, 1);
identity = repmat(struct( ...
    "surface_triangle_index", 0, ...
    "surface_triangle_nodes", zeros(1, 3), ...
    "boundary_number", 0, ...
    "boundary_name", "", ...
    "adjacent_tet_index", 0), nRows, 1);
for k = 1:nRows
    identity(k).surface_triangle_index = k;
    identity(k).surface_triangle_nodes = surfaceTriangles(k, :);
    identity(k).boundary_number = boundaryNumbers(k);
    identity(k).boundary_name = string(boundaryNames(k));
    identity(k).adjacent_tet_index = adjacentTetIndices(k);
end
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
