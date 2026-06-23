function model = assembleFirstOrderFemBemTrace(model)
%ASSEMBLEFIRSTORDERFEMBEMTRACE Assemble the P1/P1 FEM-BEM trace scaffold.

model = addFirstOrderFemBemSpaces(model);
model.topology = buildFirstOrderTopology(model);

nFem = size(model.lukas.geo.nodes, 1);
nBem = numel(model.trace.nodeIds);
traceMatrix = sparse(1:nBem, model.trace.nodeIds, 1, nBem, nFem);

model.operators = struct();
model.operators.fem = assembleLukasP1Stiffness(model, 1);
model.operators.bem = assembleGypsilabP1SurfaceMass(model);
model.operators.trace = struct( ...
    "matrix", traceMatrix, ...
    "femNodeIds", model.trace.nodeIds, ...
    "bemNodeIds", (1:nBem).', ...
    "rwgToHcurlEdgeIds", model.topology.trace.rwgToHcurlEdgeIds);
model.operators.policy = "first_order_h1_hcurl_rwg_trace_scaffold";
model.status = "operators_ready_first_order_h1_hcurl_rwg_trace";
end
