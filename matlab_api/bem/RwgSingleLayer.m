classdef RwgSingleLayer
%RWGSINGLELAYER Galerkin vector single layer on RWG dofs (static EFIE core).
%
%   L = RwgSingleLayer(rwg);
%   L.matrix    % L_ij = int int f_i(x) . f_j(y) / (4*pi*|x-y|) dS dS
%   A = L.vectorPotentialAt(c, points);   % (1/4pi) int (sum c f)(y)/r dS
%
% This is the vector-potential (partial-inductance) kernel of the static
% EFIE and of magnetostatic surface currents: for a surface current
% K = sum_d c_d f_d the field A = mu0 * vectorPotentialAt(c, x) is the
% magnetostatic vector potential (mu0 = 1 here).
%
% Assembly mirrors GalerkinSingleLayer: the test integral uses Gauss
% quadrature, the source integral is analytic on every triangle - each
% Cartesian component of an RWG function is affine in y, so
%
%   int_T f_e(y)/|x-y| dS(y) = sigma l_e/(2A) sum_k (p_k - p_opp) I1_k(x)
%
% reuses the verified single-layer P1 panel integrals with no new
% singular math.

properties (Constant)
    kernel = "vector_single_layer_f_dot_f_over_4pi_r"
    policy = "galerkin_rwg_semi_analytic_static_efie_teaching_operator"
end

properties
    rwg          % RwgSpace carrying the dofs
    quadrature   % SurfaceQuadrature for the test integral
    matrix       % dense Galerkin matrix (nDof x nDof)
end

methods
    function op = RwgSingleLayer(rwg, options)
        arguments
            rwg (1,1) RwgSpace
            options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
        end
        op.rwg = rwg;
        surface = rwg.surface;
        op.quadrature = SurfaceQuadrature(surface, options.QuadratureOrder);
        quad = op.quadrature;

        nGauss = quad.nPoints();
        nDof = rwg.ndof();
        P = {zeros(nGauss, nDof), zeros(nGauss, nDof), zeros(nGauss, nDof)};
        areas = surface.areas();
        for d = 1:nDof
            e = rwg.dofEdgeIds(d);
            for slot = 1:2
                t = rwg.edgeTriangles(e, slot);
                sigma = 3 - 2 * slot;
                Vt = surface.vtx(surface.tri(t, :), :);
                pOpp = surface.vtx(rwg.oppositeVerticesLocal(e, slot), :);
                [~, I1] = laplacePanelIntegrals(Vt, quad.points);
                contrib = sigma * rwg.dofEdgeLengths(d) / (2 * areas(t)) ...
                    * (I1 * (Vt - pOpp));
                for comp = 1:3
                    P{comp}(:, d) = P{comp}(:, d) + contrib(:, comp);
                end
            end
        end

        [Bx, By, Bz] = rwg.basisAtQuadrature(quad);
        W = spdiags(quad.weights, 0, nGauss, nGauss);
        op.matrix = (Bx.' * (W * P{1}) + By.' * (W * P{2}) + Bz.' * (W * P{3})) ...
            / (4 * pi);
    end

    function s = shape(op, dim)
        %SHAPE [nDof, nDof] of the operator (not the object array size).
        s = size(op.matrix);
        if nargin == 2
            s = s(dim);
        end
    end

    function A = vectorPotentialAt(op, coefficients, points)
        %VECTORPOTENTIALAT Analytic vector potential of an RWG current.
        arguments
            op (1,1) RwgSingleLayer
            coefficients (:,1) double
            points (:,3) double
        end
        rwg = op.rwg;
        surface = rwg.surface;
        areas = surface.areas();
        A = zeros(size(points, 1), 3);
        for d = 1:rwg.ndof()
            if coefficients(d) == 0
                continue
            end
            e = rwg.dofEdgeIds(d);
            for slot = 1:2
                t = rwg.edgeTriangles(e, slot);
                sigma = 3 - 2 * slot;
                Vt = surface.vtx(surface.tri(t, :), :);
                pOpp = surface.vtx(rwg.oppositeVerticesLocal(e, slot), :);
                [~, I1] = laplacePanelIntegrals(Vt, points);
                coef = coefficients(d) * sigma * rwg.dofEdgeLengths(d) / (2 * areas(t));
                A = A + coef * (I1 * (Vt - pOpp)) / (4 * pi);
            end
        end
    end

    function y = mtimes(op, x)
        %MTIMES op * x is the dense Galerkin apply.
        y = op.matrix * x;
    end
end
end
