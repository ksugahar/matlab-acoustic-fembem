function space = rwg(model)
%RWG First-order boundary RWG space on .vol surface triangles.

if ~isfield(model, "topology")
    model.topology = buildFirstOrderTopology(model);
end

space = struct();
space.kind = "RWG";
space.basis = "RWG0";
space.cell = "triangle-pair";
space.vtx = model.gypsilab.vtx;
space.tri = model.gypsilab.elt;
space.globalNodeIds = model.gypsilab.globalNodeIds;
space.edgesLocal = model.topology.rwg.edgesLocal;
space.edgesGlobal = model.topology.rwg.edgesGlobal;
space.dofEdgeIds = model.topology.rwg.dofEdgeIds;
space.dofEdgesLocal = model.topology.rwg.dofEdgesLocal;
space.dofEdgesGlobal = model.topology.rwg.dofEdgesGlobal;
space.edgeTriangles = model.topology.rwg.edgeTriangles;
space.oppositeVerticesLocal = model.topology.rwg.oppositeVerticesLocal;
space.hcurlEdgeIds = model.topology.rwg.hcurlEdgeIds;
end
