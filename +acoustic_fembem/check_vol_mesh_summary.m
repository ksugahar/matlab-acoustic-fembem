function check_vol_mesh_summary(volPath)
%CHECK_VOL_MESH_SUMMARY Print compact JSON metadata for a .vol mesh.

arguments
    volPath (1,1) string
end

report = acoustic_fembem.vol_mesh_summary(volPath);

summary = struct();
summary.tool = report.tool;
summary.status = report.status;
summary.ok = report.status == "ok";
summary.input = report.input;
summary.mesh_id = report.mesh_id;
summary.points = report.points;
summary.triangles = report.triangles;
summary.tets = report.tets;
summary.trace_nodes = report.trace_nodes;
summary.boundary_orientation = report.boundary_orientation;
summary.bounding_box = report.bounding_box;
summary.recommended_gui_viewer = report.recommended_gui_viewer;
summary.recommended_matlab_preview = report.recommended_matlab_preview;

disp(jsonencode(summary));
end
