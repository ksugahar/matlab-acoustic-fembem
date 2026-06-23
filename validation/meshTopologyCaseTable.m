function cases = meshTopologyCaseTable()
%MESHTOPOLOGYCASETABLE First 10 verified mesh/topology cases.

repoRoot = gypsilabRepoRoot();
fixtureRoot = fullfile(repoRoot, "fixtures", "mesh_topology");
externalCubit = "W:\00_CAE\CoreformCubit\2026_04_06_CoreformCUBIT_ExportMesh\02_Cylinder_Tet\02_Cylinder_Tet_order1.vol";

cases = [
    row("GYP-001", "unit tetra vol topology", fullfile(fixtureRoot, "unit_tetra.vol"), true)
    row("GYP-002", "four tet interior node trace", fullfile(fixtureRoot, "four_tet_interior_node.vol"), true)
    row("GYP-003", "closed tetra surface manifold", fullfile(fixtureRoot, "closed_tetra_surface_manifold.vol"), true)
    row("GYP-004", "two material tetra labels", fullfile(fixtureRoot, "two_material_tetra_labels.vol"), true)
    row("GYP-005", "boundary name propagation", fullfile(fixtureRoot, "boundary_name_propagation.vol"), true)
    row("GYP-006", "reversed surface orientation signs", fullfile(fixtureRoot, "reversed_surface_orientation.vol"), true)
    row("GYP-007", "tet edge orientation signs", fullfile(fixtureRoot, "tet_edge_orientation_signs.vol"), true)
    row("GYP-008", "quad surface rejection", fullfile(fixtureRoot, "quad_surface_rejection.vol"), false)
    row("GYP-009", "hex volume rejection", fullfile(fixtureRoot, "hex_volume_rejection.vol"), false)
    row("GYP-010", "real cubit cylinder vol intake", externalCubit, true)
];
end


function item = row(id, title, volFile, expectOk)
item = struct();
item.id = string(id);
item.title = string(title);
item.volFile = string(volFile);
item.expectOk = expectOk;
end
