function report = vol_mesh_summary(volPath)
%VOL_MESH_SUMMARY Summarize a .vol mesh for MCP and script preflight.
%
%   report = acoustic_fembem.vol_mesh_summary("unit_sphere_coarse.vol")
%
% The path may be absolute, relative to the repository, or a fixture file name
% under fixtures/mesh_topology.

arguments
    volPath (1,1) string
end

root = acoustic_fembem.repository_root();
addpath(genpath(fullfile(root, "matlab_api")));

resolved = resolveVolPath(root, volPath);
mesh = VolMesh(resolved);

bounds = struct();
bounds.min = min(mesh.vtx, [], 1);
bounds.max = max(mesh.vtx, [], 1);
bounds.span = bounds.max - bounds.min;

report = struct();
report.tool = "acoustic_fembem_vol_mesh_summary";
report.status = "ok";
report.input = volPath;
report.resolved_path = resolved;
report.mesh_id = mesh.meshId;
report.source_file_id = mesh.sourceFileId;
report.policy = mesh.policy;
report.points = size(mesh.vtx, 1);
report.triangles = size(mesh.tri, 1);
report.tets = size(mesh.tet, 1);
report.materials = mapToRows(mesh.materials, "material");
report.boundaries = mapToRows(mesh.boundaries, "boundary");
report.boundary_orientation = mesh.boundaryOrientation.boundaryOrientation;
report.trace_nodes = numel(mesh.traceNodeIds);
report.bounding_box = bounds;
report.recommended_gui_viewer = "Netgen/native .vol viewer";
report.recommended_matlab_preview = "plotVolMesh";
report.recommended_llm_preflight = "acoustic_fembem_vol_mesh_summary";
end


function path = resolveVolPath(root, requested)
requested = string(requested);
candidates = [
    requested
    string(fullfile(root, requested))
    string(fullfile(root, "fixtures", "mesh_topology", requested))
];

for k = 1:numel(candidates)
    candidate = candidates(k);
    if isfile(candidate)
        path = candidate;
        return
    end
end

error("acoustic_fembem:VolMeshNotFound", ...
    "Could not find .vol mesh: %s", requested);
end


function rows = mapToRows(map, fallbackPrefix)
keys = cell2mat(map.keys);
keys = sort(keys(:));
rows = repmat(struct("id", 0, "name", ""), numel(keys), 1);
for k = 1:numel(keys)
    id = keys(k);
    rows(k).id = id;
    if isKey(map, id)
        rows(k).name = string(map(id));
    else
        rows(k).name = string(fallbackPrefix) + "_" + id;
    end
end
end
