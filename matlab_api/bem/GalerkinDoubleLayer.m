classdef GalerkinDoubleLayer
%GALERKINDOUBLELAYER Galerkin Laplace/Helmholtz double layer on boundary P1.
%
%   K = GalerkinDoubleLayer(surface);                    % Laplace (k = 0)
%   K = GalerkinDoubleLayer(surface, "Wavenumber", k);   % Helmholtz
%   K.matrix     % K_ij = int int phi_i(x) dG_k/dn(y) (x,y) phi_j(y) dS dS
%
% The kernel is the OUTWARD-normal double layer with the e^{+ikr}
% convention,
%
%   dG_k/dn(y) = n(y).(x-y) exp(1i*k*r) (1 - 1i*k*r) / (4*pi*r^3),
%
% split like GalerkinSingleLayer: the singular Laplace part
% n.(x-y)/(4*pi*r^3) is integrated ANALYTICALLY over every source triangle
% (laplaceDoubleLayerPanelIntegrals, principal value on the diagonal - the
% +-1/2 jump belongs to the boundary integral equation, not to this
% matrix), and the smooth low-frequency-stable correction
% base * (exp(z)(1-z) - 1), z = 1i*k*r, goes through plain quadrature
% (HelmholtzKernel with SourceNormals), so the k -> 0 limit is exactly the
% Laplace operator. Stored triangle winding is corrected to the outward
% normal with the mesh orientation signs; an unknown orientation fails
% loudly.
%
% Sphere identities used by the tests (unit sphere, this kernel):
%   k = 0:  K[1] = -1/2      K[Y_l] = -1/(2*(2l+1))
%   k > 0:  K[Y_l] = 1/2 + 1i*k^2 * j_l(k) * h_l'(k)
%           (l = 0: K[1] = 1/2 - 1i*k^2 * j_0(k) * h_1(k))

properties (Constant)
    kernel = "outward_normal_double_layer_n_dot_x_minus_y_over_4pi_r3"
    policy = "galerkin_p1_semi_analytic_double_layer_principal_value"
end

properties
    surface      % SurfaceMesh carrying the boundary P1 dofs
    quadrature   % SurfaceQuadrature used for the test integral
    wavenumber   % k in the Helmholtz kernel (0 = Laplace)
    matrix       % dense Galerkin matrix (nNodes x nNodes)
end

methods
    function op = GalerkinDoubleLayer(surface, options)
        arguments
            surface (1,1) SurfaceMesh
            options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
            options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
        end
        op.surface = surface;
        op.wavenumber = options.Wavenumber;
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

        % --- smooth Helmholtz correction: plain double quadrature ---
        if op.wavenumber > 0
            parts = HelmholtzKernel(quad.points, quad.points, ...
                "Wavenumber", op.wavenumber, ...
                "SourceWeights", quad.weights, ...
                "SourceNormals", quad.outwardNormals());
            % doubleLayerSourceNormalCorrection carries weights + 1/(4*pi);
            % same-triangle pairs vanish exactly (n.(x-y) = 0 in-plane)
            op.matrix = op.matrix ...
                + Bw.' * (parts.doubleLayerSourceNormalCorrection * quad.basis);
        end
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
