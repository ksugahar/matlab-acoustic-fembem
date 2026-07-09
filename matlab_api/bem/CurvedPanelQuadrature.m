classdef CurvedPanelQuadrature
%CURVEDPANELQUADRATURE Isoparametric curved-panel Gauss quadrature (curve order 1/2/3).
%
% Superparametric surface quadrature: the panel GEOMETRY is a Lagrange curved
% triangle of curve order 1 (flat), 2 (6-node quadratic), or 3 (10-node cubic)
% -- corner nodes on the surface plus edge/interior nodes projected onto the
% true surface -- while the SOLUTION field stays P1 (values at the 3 corners).
% This is the readable "curved P1" (superparametric) element: raising the CURVE
% order (not the fes order) removes the O(h^2) flat-faceting error that caps a
% straight-panel BEM.  Because geometry error and field error are separate
% channels, on a curved boundary the geometry channel is the binding one until
% it is resolved -- so curve order is the lever (testCurvedPanelCurveOrderSweep).
%
% The single geometry knob besides curve order is a Projection function handle
% X(nx3) -> X(nx3) that snaps a point onto the surface.  Projection = @(X) X
% (the default) leaves the extra nodes at the straight lattice, every Lagrange
% map degenerates to the affine one, and this reproduces SurfaceQuadrature
% EXACTLY -- so "flat vs curved" is a one-line change and the A/B isolates
% geometry.  Curve order 1 is flat for ANY projection (only corner nodes).
%
%   proj = CurvedPanelQuadrature.sphereProjection(R);
%   cq = CurvedPanelQuadrature(surface, 7, proj);          % quadratic (default)
%   cq = CurvedPanelQuadrature(surface, 7, proj, 3);       % cubic curve order
%   fq = CurvedPanelQuadrature(surface, 7);                 % flat (== SurfaceQuadrature)
%   cq.points / cq.weights / cq.triangleIndex / cq.basis / cq.outwardNormals
%   cq.centroids()     % curved panel centroids (mapped bary [1/3 1/3 1/3])
%   cq.geomNodes       % the Lagrange geometry nodes per triangle (nTris x M x 3)

properties
    surface        % SurfaceMesh carrying the P1 nodes
    order          % Gauss points per triangle: 1, 3, or 7
    curveOrder     % geometry (Lagrange) order: 1 flat, 2 quadratic, 3 cubic
    projection     % X(nx3) -> X(nx3), snaps a point onto the surface
    geomNodes      % Lagrange geometry nodes per triangle (nTris x M x 3)
    points         % Gauss points on the curved surface ((nTris*order) x 3)
    weights        % curved quadrature weights ((nTris*order) x 1)
    triangleIndex  % source triangle of each point ((nTris*order) x 1)
    basis          % sparse (nPoints x nNodes): P1 basis at each point
end

methods
    function quad = CurvedPanelQuadrature(surface, order, projection, curveOrder)
        arguments
            surface (1,1) SurfaceMesh
            order (1,1) double {mustBeMember(order, [1 3 7])} = 7
            projection (1,1) function_handle = @(X) X
            curveOrder (1,1) double {mustBeMember(curveOrder, [1 2 3])} = 2
        end
        [bary, wfrac] = triangleGaussRule(order);
        [N, dNdu, dNdv] = CurvedPanelQuadrature.lagrangeShapes(curveOrder, bary);

        tri = surface.tri;
        vtx = surface.vtx;
        nTris = size(tri, 1);
        nNodes = size(vtx, 1);

        quad.surface = surface;
        quad.order = order;
        quad.curveOrder = curveOrder;
        quad.projection = projection;
        quad.geomNodes = CurvedPanelQuadrature.lagrangeNodes(curveOrder, vtx, tri, projection);

        quad.points = zeros(nTris * order, 3);
        quad.weights = zeros(nTris * order, 1);
        quad.triangleIndex = zeros(nTris * order, 1);

        rows = zeros(nTris * order * 3, 1);
        cols = zeros(nTris * order * 3, 1);
        vals = zeros(nTris * order * 3, 1);
        cursor = 1;
        for t = 1:nTris
            span = (t - 1) * order + (1:order);
            P = squeeze(quad.geomNodes(t, :, :));   % M x 3 geometry nodes
            Xg = N * P;                             % order x 3 curved points
            cr = cross(dNdu * P, dNdv * P, 2);      % order x 3 (2*area normal)
            detJ = sqrt(sum(cr.^2, 2));             % |x_u x x_v|
            quad.points(span, :) = Xg;
            quad.weights(span) = wfrac(:) .* 0.5 .* detJ;
            quad.triangleIndex(span) = t;
            for g = 1:order
                for k = 1:3
                    rows(cursor) = span(g);
                    cols(cursor) = tri(t, k);
                    vals(cursor) = bary(g, k);      % P1 basis == barycentric
                    cursor = cursor + 1;
                end
            end
        end
        quad.basis = sparse(rows, cols, vals, nTris * order, nNodes);
    end

    function n = nPoints(quad)
        n = numel(quad.weights);
    end

    function B = weightedBasis(quad)
        %WEIGHTEDBASIS sparse (nPoints x nNodes) with w_g * phi_i(g) entries.
        B = spdiags(quad.weights, 0, quad.nPoints(), quad.nPoints()) * quad.basis;
    end

    function c = centroids(quad)
        %CENTROIDS Curved panel centroids: the Lagrange map at bary [1/3 1/3 1/3].
        N = CurvedPanelQuadrature.lagrangeShapes(quad.curveOrder, [1 1 1] / 3);
        nTris = size(quad.geomNodes, 1);
        c = zeros(nTris, 3);
        for t = 1:nTris
            c(t, :) = N * squeeze(quad.geomNodes(t, :, :));
        end
    end

    function Nrm = outwardNormals(quad)
        %OUTWARDNORMALS Unit outward normal at each Gauss point (nPoints x 3).
        %
        % Curved geometric normal (x_u x x_v)/|x_u x x_v| corrected to outward
        % with the mesh orientation signs (fails loudly on unknown orientation).
        signs = quad.surface.orientation.triangleOrientationSignsToOutward(:);
        if any(signs == 0)
            error("CurvedPanelQuadrature:orientation", ...
                "Surface orientation is unknown for %d triangle(s).", sum(signs == 0));
        end
        [bary, ~] = triangleGaussRule(quad.order);
        [~, dNdu, dNdv] = CurvedPanelQuadrature.lagrangeShapes(quad.curveOrder, bary);
        nTris = size(quad.geomNodes, 1);
        Nrm = zeros(quad.nPoints(), 3);
        for t = 1:nTris
            span = (t - 1) * quad.order + (1:quad.order);
            P = squeeze(quad.geomNodes(t, :, :));
            cr = cross(dNdu * P, dNdv * P, 2);
            Nrm(span, :) = signs(t) * cr ./ sqrt(sum(cr.^2, 2));
        end
    end
end

methods (Static)
    function proj = sphereProjection(radius)
        %SPHEREPROJECTION Projection onto the sphere of the given radius (origin).
        arguments
            radius (1,1) double {mustBePositive} = 1.0
        end
        proj = @(X) radius * X ./ vecnorm(X, 2, 2);
    end

    function nodes = lagrangeNodes(curveOrder, vtx, tri, projection)
        %LAGRANGENODES Lagrange geometry nodes per triangle (nTris x M x 3).
        %
        % Corners stay on the surface; every higher-order lattice node is the
        % straight-triangle Lagrange position projected onto the surface.
        % Projection = @(X) X leaves them straight, so the map is exactly affine.
        b = CurvedPanelQuadrature.referenceLattice(curveOrder);   % M x 3 barycentric
        M = size(b, 1);
        nTris = size(tri, 1);
        P1 = vtx(tri(:, 1), :); P2 = vtx(tri(:, 2), :); P3 = vtx(tri(:, 3), :);
        nodes = zeros(nTris, M, 3);
        for m = 1:M
            straight = b(m, 1) * P1 + b(m, 2) * P2 + b(m, 3) * P3;
            nodes(:, m, :) = projection(straight);
        end
    end

    function b = referenceLattice(curveOrder)
        %REFERENCELATTICE Barycentric lattice of the Lagrange nodes (M x 3).
        switch curveOrder
            case 1
                b = [1 0 0; 0 1 0; 0 0 1];
            case 2
                b = [1 0 0; 0 1 0; 0 0 1
                     1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2];
            case 3
                b = [1 0 0; 0 1 0; 0 0 1
                     2/3 1/3 0; 1/3 2/3 0
                     0 2/3 1/3; 0 1/3 2/3
                     1/3 0 2/3; 2/3 0 1/3
                     1/3 1/3 1/3];
        end
    end

    function [N, dNdu, dNdv] = lagrangeShapes(curveOrder, bary)
        %LAGRANGESHAPES Lagrange triangle shapes at barycentric points (u=L2, v=L3).
        %
        % L1 = 1-u-v so dL1/du = dL1/dv = -1, dL2/du = 1, dL3/dv = 1.  Node
        % order matches referenceLattice(curveOrder).  Returns N, dN/du, dN/dv.
        L1 = bary(:, 1); L2 = bary(:, 2); L3 = bary(:, 3);
        z = zeros(size(L1)); o = ones(size(L1));
        switch curveOrder
            case 1
                N = [L1, L2, L3];
                dNdu = [-o, o, z];
                dNdv = [-o, z, o];
            case 2
                N = [L1.*(2*L1-1), L2.*(2*L2-1), L3.*(2*L3-1), ...
                     4*L1.*L2, 4*L2.*L3, 4*L3.*L1];
                dNdu = [-(4*L1-1), 4*L2-1, z, 4*(L1-L2), 4*L3, -4*L3];
                dNdv = [-(4*L1-1), z, 4*L3-1, -4*L2, 4*L2, 4*(L1-L3)];
            case 3
                % dN/du = -dN/dL1 + dN/dL2 ; dN/dv = -dN/dL1 + dN/dL3
                gL1 = [0.5*(27*L1.^2-18*L1+2), z, z, ...
                       4.5*L2.*(6*L1-1), 4.5*L2.*(3*L2-1), z, z, ...
                       4.5*L3.*(3*L3-1), 4.5*L3.*(6*L1-1), 27*L2.*L3];
                gL2 = [z, 0.5*(27*L2.^2-18*L2+2), z, ...
                       4.5*L1.*(3*L1-1), 4.5*L1.*(6*L2-1), ...
                       4.5*L3.*(6*L2-1), 4.5*L3.*(3*L3-1), z, z, 27*L1.*L3];
                gL3 = [z, z, 0.5*(27*L3.^2-18*L3+2), z, z, ...
                       4.5*L2.*(3*L2-1), 4.5*L2.*(6*L3-1), ...
                       4.5*L1.*(6*L3-1), 4.5*L1.*(3*L1-1), 27*L1.*L2];
                N = [0.5*L1.*(3*L1-1).*(3*L1-2), 0.5*L2.*(3*L2-1).*(3*L2-2), ...
                     0.5*L3.*(3*L3-1).*(3*L3-2), ...
                     4.5*L1.*L2.*(3*L1-1), 4.5*L1.*L2.*(3*L2-1), ...
                     4.5*L2.*L3.*(3*L2-1), 4.5*L2.*L3.*(3*L3-1), ...
                     4.5*L3.*L1.*(3*L3-1), 4.5*L3.*L1.*(3*L1-1), 27*L1.*L2.*L3];
                dNdu = -gL1 + gL2;
                dNdv = -gL1 + gL3;
        end
    end
end
end
