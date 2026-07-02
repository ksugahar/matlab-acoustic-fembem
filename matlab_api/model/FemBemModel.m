classdef FemBemModel
%FEMBEMMODEL Coupled first-order FEM/BEM teaching model from a .vol file.
%
%   m = FemBemModel("model.vol");
%   m.mesh        % VolMesh (vtx/tet/tri + labels + identity)
%   m.surface     % SurfaceMesh (compact boundary + row identity)
%   m.h1          % H1Space:        [K, d] = m.h1.stiffness(coef)
%   m.hcurl       % Nedelec0Space:  [M, d] = m.hcurl.mass(coef)
%   m.scalarBem   % SurfaceP1Space: [M, d] = m.scalarBem.mass()
%   m.rwg         % RwgSpace:       boundary triangle-pair dofs
%   m.trace       % TraceOperator:  g = m.trace * u
%
%   m = m.assemble();   % operators struct for the coupling manifest
%
% The constructor builds the mesh, the four first-order spaces, the trace
% operator, and the RWG-to-HCurl oriented-edge map. Matrix assembly stays in
% the named assemble method.

properties (Constant)
    spacesPolicy = "first_order_h1_p1_hcurl_nedelec0_bem_p1_rwg_only"
    operatorsPolicy = "first_order_h1_hcurl_rwg_trace_scaffold"
end

properties
    mesh                % VolMesh
    surface             % SurfaceMesh
    h1                  % H1Space (volume P1)
    hcurl               % Nedelec0Space (volume edge)
    scalarBem           % SurfaceP1Space (boundary P1)
    rwg                 % RwgSpace (boundary triangle-pair)
    trace               % TraceOperator (H1 -> boundary P1)
    rwgToHcurlEdgeIds   % volume edge id of every RWG dof edge
    operators           % [] until assemble(): fem/bem/trace operator bundle
    status              % readable model state
end

methods
    function model = FemBemModel(volFile)
        arguments
            volFile (1,1) string
        end
        model.mesh = VolMesh(volFile);
        model.surface = model.mesh.boundary();
        model.h1 = H1Space(model.mesh);
        model.hcurl = Nedelec0Space(model.mesh);
        model.scalarBem = SurfaceP1Space(model.surface);
        model.rwg = RwgSpace(model.surface);
        model.trace = TraceOperator(model.mesh);
        model.rwgToHcurlEdgeIds = model.rwg.hcurlEdgeIds(model.hcurl);
        model.operators = [];
        model.status = "vol_ready_first_order_h1_hcurl_rwg";
    end

    function model = assemble(model)
        %ASSEMBLE Assemble the P1/P1 FEM-BEM operator scaffold.
        %
        % fem:   volume P1 stiffness with unit coefficient
        % bem:   boundary P1 surface mass
        % trace: the one-hot TraceOperator plus the RWG-to-HCurl edge map
        [~, femDetail] = model.h1.stiffness(1);
        [~, bemDetail] = model.scalarBem.mass();

        model.operators = struct();
        model.operators.fem = femDetail;
        model.operators.bem = bemDetail;
        model.operators.trace = model.trace;
        model.operators.rwgToHcurlEdgeIds = model.rwgToHcurlEdgeIds;
        model.operators.policy = model.operatorsPolicy;
        model.status = "operators_ready_first_order_h1_hcurl_rwg_trace";
    end

    function catalog = spaceCatalog(model)
        %SPACECATALOG The fixed first-order space family of this lane.
        catalog = struct();
        catalog.h1 = spaceRow(model.h1);
        catalog.hcurl = spaceRow(model.hcurl);
        catalog.scalarBem = spaceRow(model.scalarBem);
        catalog.rwg = spaceRow(model.rwg);
        catalog.policy = model.spacesPolicy;
    end
end
end


function row = spaceRow(space)
%SPACEROW Family/order/cell/kind metadata of one space object.

row = struct( ...
    "family", space.family, ...
    "order", space.order, ...
    "cell", space.cell, ...
    "kind", space.basis);
end
