classdef GalerkinDoubleLayer
%GALERKINDOUBLELAYER Galerkin Laplace double-layer operator on boundary P1.
%
%   K = GalerkinDoubleLayer(surface);
%   K.matrix     % K_ij = int int phi_i(x) dG/dn(y) (x,y) phi_j(y) dS dS
%
% The kernel is the OUTWARD-normal double layer
%
%   dG/dn(y) = n(y) . (x - y) / (4*pi*|x - y|^3),
%
% assembled like GalerkinSingleLayer: test integral by Gauss quadrature,
% source integral analytic on every triangle
% (laplaceDoubleLayerPanelIntegrals, principal value on the diagonal - the
% +-1/2 jump belongs to the boundary integral equation, not to this
% matrix). Stored triangle winding is corrected to the outward normal with
% the mesh orientation signs; an unknown orientation fails loudly.
%
% Sphere identities used by the tests (unit sphere, this kernel):
%   K[1]   = -1/2        K[Y_1] = -1/6      (K[Y_l] = -1/(2*(2l+1)))

properties (Constant)
    kernel = "outward_normal_double_layer_n_dot_x_minus_y_over_4pi_r3"
    policy = "galerkin_p1_semi_analytic_double_layer_principal_value"
end

properties
    surface      % SurfaceMesh carrying the boundary P1 dofs
    quadrature   % SurfaceQuadrature used for the test integral
    matrix       % dense Galerkin matrix (nNodes x nNodes)
end

methods
    function op = GalerkinDoubleLayer(surface, options)
        arguments
            surface (1,1) SurfaceMesh
            options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
        end
        op.surface = surface;
        op.quadrature = SurfaceQuadrature(surface, options.QuadratureOrder);
        quad = op.quadrature;

        signs = surface.orientation.triangleOrientationSignsToOutward(:);
        if any(signs == 0)
            error("GalerkinDoubleLayer:orientation", ...
                "Surface orientation is unknown for %d triangle(s); cannot fix the outward normal.", ...
                sum(signs == 0));
        end

        nGauss = quad.nPoints();
        nNodes = size(surface.vtx, 1);
        P = zeros(nGauss, nNodes);
        tri = surface.tri;
        vtx = surface.vtx;
        for t = 1:size(tri, 1)
            [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t, :), :), quad.points);
            P(:, tri(t, :)) = P(:, tri(t, :)) + signs(t) * J1;
        end
        Bw = quad.weightedBasis();
        op.matrix = Bw.' * P / (4 * pi);
    end

    function s = shape(op, dim)
        %SHAPE [nNodes, nNodes] of the operator (not the object array size).
        s = size(op.matrix);
        if nargin == 2
            s = s(dim);
        end
    end

    function p = apply(op, g)
        %APPLY Galerkin double-layer action on a P1 trace.
        p = op.matrix * g;
    end

    function p = mtimes(op, g)
        %MTIMES op * g is the same dense apply.
        p = op.matrix * g;
    end
end
end
