function space = h1(model)
%H1 First-order H1 space on .vol tetrahedra.

space = struct();
space.kind = "H1";
space.basis = "P1";
space.cell = "tetrahedron";
space.nodes = model.lukas.geo.nodes;
space.tet = model.lukas.geo.conn_matrix;
space.traceNodeIds = model.trace.nodeIds;
space.stiffness = @assemble;

    function out = assemble(materialCoef)
        if nargin == 0
            materialCoef = 1;
        end
        out = assembleLukasP1Stiffness(model, materialCoef);
    end
end
