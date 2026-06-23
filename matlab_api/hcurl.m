function space = hcurl(model)
%HCURL First-order HCurl space on .vol tetrahedra.

if ~isfield(model, "topology")
    model.topology = buildFirstOrderTopology(model);
end

space = struct();
space.kind = "HCurl";
space.basis = "Nedelec0";
space.cell = "tetrahedron";
space.edges = model.topology.hcurl.edges;
space.tetEdges = model.topology.hcurl.tetEdges;
space.tetEdgeSigns = model.topology.hcurl.tetEdgeSigns;
space.traceRwgEdgeIds = model.topology.trace.rwgToHcurlEdgeIds;
end
