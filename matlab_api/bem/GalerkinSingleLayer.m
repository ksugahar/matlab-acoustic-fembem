classdef GalerkinSingleLayer
%GALERKINSINGLELAYER Galerkin single-layer operator on boundary P1.
%
%   V = GalerkinSingleLayer(surface);                    % Laplace (k = 0)
%   V = GalerkinSingleLayer(surface, "Wavenumber", k);   % Helmholtz
%   V.matrix          % (nNodes x nNodes), V_ij = int int phi_i G_k phi_j
%   p = V * q;
%
% The Galerkin double integral
%
%   V_ij = int_S int_S phi_i(x) G_k(x,y) phi_j(y) dS(y) dS(x),
%   G_k  = 1/(4*pi*r) + (exp(1i*k*r) - 1)/(4*pi*r)
%
% is split exactly like Gypsilab's integral + regularize pair:
%
%   test integral  : Gauss quadrature (SurfaceQuadrature, 1/3/7 points)
%   singular part  : the Laplace kernel is integrated ANALYTICALLY over
%                    every source triangle (laplacePanelIntegrals) at each
%                    test Gauss point - the semi-analytic correction
%                    applied to ALL pairs, not only near ones
%   smooth part    : the low-frequency-stable Helmholtz correction is
%                    regular, so plain quadrature (HelmholtzKernel) is fine
%
% For k = 0 the operator is purely real and the smooth part vanishes; the
% k -> 0 limit is therefore exact by construction.

properties (Constant)
    kernel = "exp(1i*k*r)/(4*pi*r)"
    policy = "galerkin_p1_semi_analytic_laplace_plus_quadrature_correction"
end

properties
    surface      % SurfaceMesh carrying the boundary P1 dofs
    quadrature   % SurfaceQuadrature used for test and smooth integrals
    wavenumber   % k in the Helmholtz kernel (0 = Laplace)
    matrix       % dense Galerkin matrix (nNodes x nNodes)
end

methods
    function op = GalerkinSingleLayer(surface, options)
        arguments
            surface (1,1) SurfaceMesh
            options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
            options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
        end
        op.surface = surface;
        op.wavenumber = options.Wavenumber;
        op.quadrature = SurfaceQuadrature(surface, options.QuadratureOrder);
        quad = op.quadrature;

        % --- singular Laplace part: analytic source integral per triangle ---
        % P(g, i) accumulates int_T lambda_i(y)/|x_g - y| dS(y) over all
        % source triangles T that carry dof i.
        nGauss = quad.nPoints();
        nNodes = size(surface.vtx, 1);
        P = zeros(nGauss, nNodes);
        tri = surface.tri;
        vtx = surface.vtx;
        for t = 1:size(tri, 1)
            [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), quad.points);
            P(:, tri(t, :)) = P(:, tri(t, :)) + I1;
        end
        Bw = quad.weightedBasis();
        op.matrix = Bw.' * P / (4 * pi);

        % --- smooth Helmholtz correction: plain double quadrature ---
        if op.wavenumber > 0
            parts = HelmholtzKernel(quad.points, quad.points, ...
                "Wavenumber", op.wavenumber, ...
                "SourceWeights", quad.weights);
            % singleLayerCorrection already carries source weights + 1/(4*pi)
            op.matrix = op.matrix + Bw.' * (parts.singleLayerCorrection * quad.basis);
        end
    end

    function s = shape(op, dim)
        %SHAPE [nNodes, nNodes] of the operator (not the object array size).
        s = size(op.matrix);
        if nargin == 2
            s = s(dim);
        end
    end

    function p = apply(op, q)
        %APPLY Galerkin single-layer action on a P1 density.
        p = op.matrix * q;
    end

    function p = mtimes(op, q)
        %MTIMES op * q is the same dense apply.
        p = op.matrix * q;
    end
end
end
