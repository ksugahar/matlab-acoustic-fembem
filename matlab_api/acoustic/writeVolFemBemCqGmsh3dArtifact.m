function artifact = writeVolFemBemCqGmsh3dArtifact(volFile, options)
%WRITEVOLFEMBEMCQGMSH3DARTIFACT Export .vol P1 FEM/BEM CQ as Gmsh v4.1.
%
%   artifact = writeVolFemBemCqGmsh3dArtifact("mesh.vol") reads a Netgen
%   first-order tri/tet .vol mesh, runs the readable P1 volume-FEM / P1
%   boundary-BEM Johnson-Nedelec convolution-quadrature solver, and writes a
%   Gmsh MSH v4.1 tetrahedral volume carrying time-stepped NodeData:
%
%       interior_pressure      on all P1 volume nodes
%       boundary_density       on boundary nodes, zero on interior nodes
%
%   The writer deliberately keeps the API small and Gypsilab-like: mesh in,
%   operators implied by the model, artifact out.  The high-order impedance
%   boundary policy is recorded as the open-boundary lane metadata; the
%   numerical radiation condition in this function is the coupled BEM row.

arguments
    volFile (1,1) string = ""
    options.OutputBase (1,1) string = ""
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 3)} = 16
    options.TimeStep (1,1) double {mustBePositive} = 0.04
    options.ExteriorSoundSpeed (1,1) double {mustBePositive} = 1.0
    options.InteriorSoundSpeed (1,1) double {mustBePositive} = 1.0
    options.Method (1,1) string {mustBeMember(options.Method, ["BDF1", "BDF2"])} = "BDF1"
    options.CouplingForm (1,1) string {mustBeMember(options.CouplingForm, ...
        ["JohnsonNedelec", "SingleLayerTeaching"])} = "JohnsonNedelec"
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 1
    options.HighOrderImpedanceBoundary (1,1) logical = true
    options.BoundaryPolicy (1,1) string = "bem_radiation_closure_with_high_order_impedance_boundary_lane"
    options.ExteriorCutPlane (1,1) logical = false
    options.ExteriorCutPlaneGrid (1,1) double {mustBeInteger, mustBeGreaterThan(options.ExteriorCutPlaneGrid, 3)} = 41
    options.ExteriorCutPlaneRadiusFactor (1,1) double {mustBeGreaterThan(options.ExteriorCutPlaneRadiusFactor, 1)} = 2.8
    options.ExteriorCutPlaneInnerRadiusFactor (1,1) double {mustBeGreaterThan(options.ExteriorCutPlaneInnerRadiusFactor, 1)} = 1.02
end

if strlength(volFile) == 0
    volFile = defaultFixture();
end
if strlength(options.OutputBase) == 0
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    options.OutputBase = fullfile(tempdir, "vol_fembem_cq_gmsh3d_" + stamp);
end

totalTimer = tic;
model = FemBemModel(volFile);
mesh = model.mesh;
surface = model.surface;
exteriorViz = makeExteriorCutPlane(mesh.vtx, options);
solverArgs = { ...
    "NumTime", options.NumTime, ...
    "TimeStep", options.TimeStep, ...
    "ExteriorSoundSpeed", options.ExteriorSoundSpeed, ...
    "InteriorSoundSpeed", options.InteriorSoundSpeed, ...
    "Method", options.Method, ...
    "CouplingForm", options.CouplingForm, ...
    "QuadratureOrder", options.QuadratureOrder};
if exteriorViz.enabled
    solverArgs = [solverArgs, {"ObservationPoints", exteriorViz.points}];
end
solveTimer = tic;
result = volFemBemCoupledConvolutionQuadrature(volFile, solverArgs{:});
solveSeconds = toc(solveTimer);
exteriorViz.initialTimeStep = 0;
if exteriorViz.enabled
    exteriorViz.initialTimeStep = max(0, min(numel(result.time) - 1, round(0.45 * numel(result.time))));
    exteriorViz.displayRange = max(eps, 0.6 * result.summary.max_abs_exterior_pressure);
end

outBase = options.OutputBase;
mshPath = outBase + ".msh";
geoPath = outBase + ".geo";
optPath = outBase + ".opt";
geoOptPath = geoPath + ".opt";
mshOptPath = mshPath + ".opt";
matPath = outBase + ".mat";
buildJsonPath = outBase + ".build.json";

writeTimer = tic;
writeGmsh41(mshPath, mesh, surface, result, exteriorViz);
writeSeconds = toc(writeTimer);

geoTimer = tic;
writeGmshCompanion(geoPath, mshPath, exteriorViz);
geoSeconds = toc(geoTimer);

optTimer = tic;
writeGmshOptions(optPath, exteriorViz);
copyfile(optPath, geoOptPath, "f");
writeGmshMeshOptions(mshOptPath);
optSeconds = toc(optTimer);

matTimer = tic;
save(matPath, "result", "mesh", "surface", "-v7.3");
matSeconds = toc(matTimer);

jsonTimer = tic;
artifact = struct();
artifact.schema = "matlab-acoustic-fembem.vol-fembem-cq-gmsh3d-build.v1";
artifact.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
artifact.matlab_version = string(version);
artifact.vol_file = string(volFile);
artifact.mesh_id = result.meshId;
artifact.mesh_source_id = result.meshSourceId;
artifact.gmsh_msh = mshPath;
artifact.gmsh_geo = geoPath;
artifact.gmsh_opt = optPath;
artifact.gmsh_geo_opt = geoOptPath;
artifact.gmsh_msh_opt = mshOptPath;
artifact.mat_result = matPath;
artifact.gmsh_msh_version = "4.1";
artifact.data_kind = "3D tetrahedral P1 FEM pressure plus P1 boundary BEM density NodeData";
if exteriorViz.enabled
    artifact.data_kind = artifact.data_kind + " plus exterior acoustic cut-plane propagation NodeData";
end
artifact.method = result.method;
artifact.coupling_form = result.couplingForm;
artifact.boundary_policy = options.BoundaryPolicy;
artifact.high_order_impedance_boundary_lane = options.HighOrderImpedanceBoundary;
artifact.boundary_closure_note = "Numerical radiation is the coupled P1 BEM row; the high-order impedance boundary is recorded as the Radia open-boundary policy lane, not Kelvin.";
artifact.nodes = size(mesh.vtx, 1);
artifact.tetrahedra = size(mesh.tet, 1);
artifact.boundary_nodes = size(surface.vtx, 1);
artifact.boundary_triangles = size(surface.tri, 1);
artifact.exterior_cut_plane_enabled = exteriorViz.enabled;
artifact.exterior_cut_plane_nodes = size(exteriorViz.points, 1);
artifact.exterior_cut_plane_triangles = size(exteriorViz.tri, 1);
artifact.exterior_cut_plane_y = exteriorViz.y;
artifact.exterior_cut_plane_radius = exteriorViz.radius;
artifact.frames = numel(result.time);
artifact.time_start = result.time(1);
artifact.time_end = result.time(end);
artifact.max_abs_interior_pressure = result.summary.max_abs_interior_pressure;
artifact.max_abs_boundary_density = result.summary.max_abs_boundary_density;
artifact.max_abs_exterior_pressure = result.summary.max_abs_exterior_pressure;
artifact.max_relative_residual = result.summary.max_relative_residual;
artifact.max_condition_number = result.summary.max_condition_number;
artifact.max_double_layer_frobenius_norm = result.summary.max_double_layer_frobenius_norm;
artifact.checks = result.checks;
artifact.checks.high_order_impedance_boundary_lane_recorded = options.HighOrderImpedanceBoundary;
artifact.checks.not_kelvin_boundary = ~contains(lower(options.BoundaryPolicy), "kelvin");
artifact.checks.gmsh_v41 = true;
artifact.checks.tetrahedral_volume_mesh = size(mesh.tet, 1) > 0;
artifact.checks.node_data_pressure_written = true;
artifact.checks.node_data_boundary_density_written = true;
artifact.checks.exterior_cut_plane_propagation_view_written = ...
    ~options.ExteriorCutPlane || exteriorViz.enabled;
artifact.checks.gmsh_geo_opt_exact_autoload_written = isfile(geoOptPath);
artifact.checks.gmsh_msh_opt_exact_autoload_written = isfile(mshOptPath);
artifact.timing = struct( ...
    "coupled_cq_solve_seconds", solveSeconds, ...
    "gmsh_write_seconds", writeSeconds, ...
    "geo_write_seconds", geoSeconds, ...
    "opt_write_seconds", optSeconds, ...
    "mat_write_seconds", matSeconds, ...
    "json_write_seconds", toc(jsonTimer), ...
    "total_seconds", toc(totalTimer));
jsonText = jsonencode(artifact);
writeString(buildJsonPath, jsonText);
end


function viz = makeExteriorCutPlane(nodes, options)
viz = struct( ...
    "enabled", false, ...
    "points", zeros(0, 3), ...
    "tri", zeros(0, 3), ...
    "y", NaN, ...
    "radius", NaN, ...
    "grid", options.ExteriorCutPlaneGrid, ...
    "displayRange", 1.0, ...
    "initialTimeStep", 0);
if ~options.ExteriorCutPlane
    return
end

center = mean(nodes, 1);
radius = max(vecnorm(nodes - center, 2, 2));
if ~(isfinite(radius) && radius > 0)
    error("writeVolFemBemCqGmsh3dArtifact:exteriorCutPlane", ...
        "Cannot build exterior cut-plane for a degenerate mesh.");
end

n = options.ExteriorCutPlaneGrid;
span = options.ExteriorCutPlaneRadiusFactor * radius;
x = linspace(center(1) - span, center(1) + span, n);
z = linspace(center(3) - span, center(3) + span, n);
[X, Z] = meshgrid(x, z);
Y = center(2) * ones(size(X));
allPoints = [X(:), Y(:), Z(:)];

outside = vecnorm(allPoints - center, 2, 2) >= ...
    options.ExteriorCutPlaneInnerRadiusFactor * radius;
gridIndex = reshape(1:size(allPoints, 1), n, n);
tri = zeros(0, 3);
for iz = 1:n-1
    for ix = 1:n-1
        ids = [gridIndex(iz, ix), gridIndex(iz, ix + 1), ...
            gridIndex(iz + 1, ix + 1), gridIndex(iz + 1, ix)];
        if all(outside(ids))
            tri(end + 1, :) = ids([1 2 3]); %#ok<AGROW>
            tri(end + 1, :) = ids([1 3 4]); %#ok<AGROW>
        end
    end
end
if isempty(tri)
    error("writeVolFemBemCqGmsh3dArtifact:exteriorCutPlane", ...
        "Exterior cut-plane grid did not leave any exterior triangles.");
end

used = unique(tri(:));
remap = zeros(size(allPoints, 1), 1);
remap(used) = 1:numel(used);
viz.enabled = true;
viz.points = allPoints(used, :);
viz.tri = remap(tri);
viz.y = center(2);
viz.radius = radius;
end


function writeGmsh41(path, mesh, surface, result, exteriorViz)
nodes = mesh.vtx;
tets = orientTetsPositive(nodes, mesh.tet);
tris = orientSurfaceTriangles(surface.triGlobal);
nVolumeNodes = size(nodes, 1);
nExteriorNodes = size(exteriorViz.points, 1);
totalNodes = nVolumeNodes + nExteriorNodes;
boundaryValues = zeros(totalNodes, numel(result.time));
boundaryValues(surface.volNodeIds, :) = result.boundaryDensity.';
interiorValues = zeros(totalNodes, numel(result.time));
interiorValues(1:nVolumeNodes, :) = result.interiorPressure.';
if exteriorViz.enabled
    exteriorValues = zeros(totalNodes, numel(result.time));
    exteriorValues(nVolumeNodes + (1:nExteriorNodes), :) = result.exteriorPressure.';
end

allNodes = nodes;
if exteriorViz.enabled
    allNodes = [allNodes; exteriorViz.points]; %#ok<AGROW>
end
xmin = min(allNodes(:, 1)); xmax = max(allNodes(:, 1));
ymin = min(allNodes(:, 2)); ymax = max(allNodes(:, 2));
zmin = min(allNodes(:, 3)); zmax = max(allNodes(:, 3));
volXmin = min(nodes(:, 1)); volXmax = max(nodes(:, 1));
volYmin = min(nodes(:, 2)); volYmax = max(nodes(:, 2));
volZmin = min(nodes(:, 3)); volZmax = max(nodes(:, 3));

fid = fopen(path, "w");
if fid < 0
    error("writeVolFemBemCqGmsh3dArtifact:file", ...
        "Could not open Gmsh output: %s", path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n");
numPhysical = 2 + double(exteriorViz.enabled);
fprintf(fid, "$PhysicalNames\n%d\n", numPhysical);
fprintf(fid, "2 1 ""boundary_p1_bem_surface""\n");
fprintf(fid, "3 2 ""vol_p1_fem_acoustic_domain""\n");
if exteriorViz.enabled
    fprintf(fid, "2 3 ""exterior_pressure_y0_cutplane""\n");
end
fprintf(fid, "$EndPhysicalNames\n");
numSurfaces = 1 + double(exteriorViz.enabled);
fprintf(fid, "$Entities\n0 0 %d 1\n", numSurfaces);
fprintf(fid, "1 %.17g %.17g %.17g %.17g %.17g %.17g 1 1 0\n", ...
    volXmin, volYmin, volZmin, volXmax, volYmax, volZmax);
if exteriorViz.enabled
    fprintf(fid, "2 %.17g %.17g %.17g %.17g %.17g %.17g 1 3 0\n", ...
        xmin, exteriorViz.y, zmin, xmax, exteriorViz.y, zmax);
end
fprintf(fid, "1 %.17g %.17g %.17g %.17g %.17g %.17g 1 2 0\n", ...
    volXmin, volYmin, volZmin, volXmax, volYmax, volZmax);
fprintf(fid, "$EndEntities\n");

numNodeBlocks = 1 + double(exteriorViz.enabled);
fprintf(fid, "$Nodes\n%d %d 1 %d\n", numNodeBlocks, totalNodes, totalNodes);
fprintf(fid, "3 1 0 %d\n", nVolumeNodes);
for node = 1:nVolumeNodes
    fprintf(fid, "%d\n", node);
end
for node = 1:nVolumeNodes
    fprintf(fid, "%.17g %.17g %.17g\n", nodes(node, 1), nodes(node, 2), nodes(node, 3));
end
if exteriorViz.enabled
    fprintf(fid, "2 2 0 %d\n", nExteriorNodes);
    for node = 1:nExteriorNodes
        fprintf(fid, "%d\n", nVolumeNodes + node);
    end
    for node = 1:nExteriorNodes
        fprintf(fid, "%.17g %.17g %.17g\n", ...
            exteriorViz.points(node, 1), exteriorViz.points(node, 2), exteriorViz.points(node, 3));
    end
end
fprintf(fid, "$EndNodes\n");

planeTris = exteriorViz.tri + nVolumeNodes;
numElements = size(tris, 1) + size(planeTris, 1) + size(tets, 1);
numElementBlocks = 2 + double(exteriorViz.enabled);
fprintf(fid, "$Elements\n%d %d 1 %d\n", numElementBlocks, numElements, numElements);
fprintf(fid, "2 1 2 %d\n", size(tris, 1));
for elem = 1:size(tris, 1)
    fprintf(fid, "%d %d %d %d\n", elem, tris(elem, 1), tris(elem, 2), tris(elem, 3));
end
nextTag = size(tris, 1) + 1;
if exteriorViz.enabled
    fprintf(fid, "2 2 2 %d\n", size(planeTris, 1));
    for elem = 1:size(planeTris, 1)
        tag = nextTag + elem - 1;
        fprintf(fid, "%d %d %d %d\n", tag, ...
            planeTris(elem, 1), planeTris(elem, 2), planeTris(elem, 3));
    end
    nextTag = nextTag + size(planeTris, 1);
end
fprintf(fid, "3 1 4 %d\n", size(tets, 1));
for elem = 1:size(tets, 1)
    tag = nextTag + elem - 1;
    fprintf(fid, "%d %d %d %d %d\n", tag, tets(elem, 1), tets(elem, 2), tets(elem, 3), tets(elem, 4));
end
fprintf(fid, "$EndElements\n");

if exteriorViz.enabled
    writeNodeData(fid, "exterior_pressure_y0_cutplane", result.time, exteriorValues);
end
writeNodeData(fid, "interior_pressure", result.time, interiorValues);
writeNodeData(fid, "boundary_density", result.time, boundaryValues);
clear cleanup
end


function writeNodeData(fid, name, time, values)
for step = 1:numel(time)
    fprintf(fid, "$NodeData\n");
    fprintf(fid, "1\n");
    fprintf(fid, """%s""\n", name);
    fprintf(fid, "1\n");
    fprintf(fid, "%.17g\n", time(step));
    fprintf(fid, "3\n");
    fprintf(fid, "%d\n", step - 1);
    fprintf(fid, "1\n");
    fprintf(fid, "%d\n", size(values, 1));
    for node = 1:size(values, 1)
        fprintf(fid, "%d %.17g\n", node, values(node, step));
    end
    fprintf(fid, "$EndNodeData\n");
end
end


function writeGmshCompanion(geoPath, mshPath, exteriorViz)
[~, mshBase, mshExt] = fileparts(mshPath);
fid = fopen(geoPath, "w");
if fid < 0
    error("writeVolFemBemCqGmsh3dArtifact:file", ...
        "Could not open Gmsh companion: %s", geoPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Gmsh companion for .vol P1 volume-FEM / P1 boundary-BEM CQ output.\n");
fprintf(fid, "Merge ""%s%s"";\n", mshBase, mshExt);
fprintf(fid, "General.InitialModule = 5;\n");
writeGmshPositiveYCamera(fid);
writeGmshDisplayOptions(fid, exteriorViz);
clear cleanup
end


function writeGmshOptions(optPath, exteriorViz)
fid = fopen(optPath, "w", "n", "UTF-8");
if fid < 0
    error("writeVolFemBemCqGmsh3dArtifact:file", ...
        "Could not open Gmsh options file: %s", optPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Gmsh display options for .vol P1 FEM/BEM CQ output.\n");
fprintf(fid, "General.InitialModule = 5;\n");
writeGmshPositiveYCamera(fid);
writeGmshDisplayOptions(fid, exteriorViz);
clear cleanup
end


function writeGmshMeshOptions(optPath)
fid = fopen(optPath, "w", "n", "UTF-8");
if fid < 0
    error("writeVolFemBemCqGmsh3dArtifact:file", ...
        "Could not open Gmsh mesh options file: %s", optPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "// Gmsh mesh-inspection options for .msh raw mesh/data.\n");
fprintf(fid, "General.InitialModule = 1;\n");
writeGmshPositiveYCamera(fid);
fprintf(fid, "General.Color.Background = {255,255,255};\n");
fprintf(fid, "General.Color.Foreground = {0,0,0};\n");
fprintf(fid, "General.Color.Text = {0,0,0};\n");
fprintf(fid, "Mesh.SurfaceFaces = 1;\n");
fprintf(fid, "Mesh.SurfaceEdges = 1;\n");
fprintf(fid, "Mesh.VolumeEdges = 0;\n");
fprintf(fid, "Mesh.ColorCarousel = 2;\n");
fprintf(fid, "Mesh.Light = 1;\n");
fprintf(fid, "Mesh.LightTwoSide = 1;\n");
fprintf(fid, "Mesh.Clip = 0;\n");
fprintf(fid, "View[0].Visible = 0;\n");
fprintf(fid, "View[1].Visible = 0;\n");
fprintf(fid, "View[2].Visible = 0;\n");
clear cleanup
end


function writeGmshDisplayOptions(fid, exteriorViz)
fprintf(fid, "General.Color.Background = {255,255,255};\n");
fprintf(fid, "General.Color.Foreground = {0,0,0};\n");
fprintf(fid, "General.Color.Text = {0,0,0};\n");
fprintf(fid, "Mesh.SurfaceFaces = 0;\n");
fprintf(fid, "Mesh.SurfaceEdges = 1;\n");
fprintf(fid, "Mesh.VolumeEdges = 0;\n");
fprintf(fid, "Mesh.Light = 1;\n");
fprintf(fid, "Mesh.LightTwoSide = 1;\n");
fprintf(fid, "Mesh.Clip = 0;\n");
fprintf(fid, "General.Clip0A = 0;\n");
fprintf(fid, "General.Clip0B = -1;\n");
fprintf(fid, "General.Clip0C = 0;\n");
fprintf(fid, "General.Clip0D = 0;\n");
fprintf(fid, "General.ClipFactor = 5;\n");
fprintf(fid, "General.ClipWholeElements = 0;\n");
fprintf(fid, "General.ClipOnlyVolume = 1;\n");
fprintf(fid, "General.ClipOnlyDrawIntersectingVolume = 0;\n");
fprintf(fid, "View[0].Visible = 1;\n");
fprintf(fid, "View[0].Clip = %d;\n", double(~exteriorViz.enabled));
fprintf(fid, "View[0].TimeStep = %d;\n", exteriorViz.initialTimeStep);
fprintf(fid, "View[0].ShowElement = %d;\n", double(exteriorViz.enabled));
fprintf(fid, "View[0].ShowScale = 1;\n");
fprintf(fid, "View[0].IntervalsType = 2;\n");
fprintf(fid, "View[0].NbIso = 32;\n");
fprintf(fid, "View[0].ColormapNumber = 2;\n");
fprintf(fid, "View[0].Light = 0;\n");
if exteriorViz.enabled
    fprintf(fid, "View[0].RangeType = 2;\n");
    fprintf(fid, "View[0].CustomMin = %.17g;\n", -exteriorViz.displayRange);
    fprintf(fid, "View[0].CustomMax = %.17g;\n", exteriorViz.displayRange);
end
fprintf(fid, "View[1].Visible = 0;\n");
fprintf(fid, "View[1].ShowElement = 1;\n");
fprintf(fid, "View[1].IntervalsType = 2;\n");
fprintf(fid, "View[2].Visible = 0;\n");
fprintf(fid, "PostProcessing.AnimationDelay = 0.08;\n");
fprintf(fid, "PostProcessing.AnimationCycle = 1;\n");
fprintf(fid, "PostProcessing.AnimationStep = 1;\n");
end


function writeGmshPositiveYCamera(fid)
fprintf(fid, "General.Trackball = 0;\n");
fprintf(fid, "General.RotationX = 75;\n");
fprintf(fid, "General.RotationY = 0;\n");
fprintf(fid, "General.RotationZ = -20;\n");
fprintf(fid, "General.RotationCenterGravity = 1;\n");
fprintf(fid, "General.ScaleX = 1.15;\n");
fprintf(fid, "General.ScaleY = 1.15;\n");
fprintf(fid, "General.ScaleZ = 1.15;\n");
fprintf(fid, "General.TranslationX = 0;\n");
fprintf(fid, "General.TranslationY = 0;\n");
fprintf(fid, "General.TranslationZ = 0;\n");
end


function tri = orientSurfaceTriangles(tri)
%ORIENTSURFACETRIANGLES Keep the stored .vol boundary triangles as the BEM surface.
%
% The .vol reader already records orientation diagnostics; the Gmsh artifact
% only needs the visual/BEM surface connectivity.
end


function tet = orientTetsPositive(nodes, tet)
for e = 1:size(tet, 1)
    X = nodes(tet(e, :), :);
    signedSixVolume = det([X(2, :) - X(1, :); X(3, :) - X(1, :); X(4, :) - X(1, :)]);
    if signedSixVolume < 0
        tet(e, [3 4]) = tet(e, [4 3]);
    end
end
end


function writeString(path, text)
fid = fopen(path, "w", "n", "UTF-8");
if fid < 0
    error("writeVolFemBemCqGmsh3dArtifact:file", "Could not open %s", path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", char(text));
clear cleanup
end


function volFile = defaultFixture()
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", ...
    "four_tet_interior_node.vol"));
end
