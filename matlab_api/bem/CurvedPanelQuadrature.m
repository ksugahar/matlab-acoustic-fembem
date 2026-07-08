classdef CurvedPanelQuadrature
%CURVEDPANELQUADRATURE Isoparametric (quadratic) curved-panel Gauss quadrature.
%
% Superparametric surface quadrature: the panel GEOMETRY is a curved 6-node
% (quadratic) triangle -- the 3 corner nodes on the surface plus 3 edge
% nodes obtained by projecting each straight edge midpoint onto the true
% surface -- while the SOLUTION field stays P1 (values at the 3 corners).
% This is the readable "curved P1" (superparametric) element: curving the
% geometry alone removes the O(h^2) flat-faceting error that caps a
% straight-panel BEM (testHelmholtzScattering documents that its analytic
% deviations are "faceted-geometry dominated").
%
% The single knob is a Projection function handle X(nx3) -> X(nx3) that
% snaps a point onto the surface.  The default Projection = @(X) X leaves
% the edge nodes at the straight midpoints, the quadratic map degenerates
% to the affine one, and this reproduces SurfaceQuadrature EXACTLY -- so
% "flat vs curved" is a one-line change and the A/B isolates geometry.
%
%   proj = CurvedPanelQuadrature.sphereProjection(R);
%   cq = CurvedPanelQuadrature(surface, 7, proj);   % curved isoparametric
%   fq = CurvedPanelQuadrature(surface, 7);          % flat (== SurfaceQuadrature)
%   cq.points          % Gauss points on the curved surface ((nTris*order) x 3)
%   cq.weights         % w * 0.5 * |x_u x x_v|  (curved surface Jacobian)
%   cq.basis           % sparse (nPoints x nNodes) P1 values at the points
%   cq.triangleIndex   % source triangle of each point
%   cq.centroids()     % curved panel centroids (mapped bary [1/3 1/3 1/3])
%   cq.outwardNormals()% unit outward normal at each Gauss point
%   cq.nodes6          % the 6 geometry nodes per triangle (nTris x 6 x 3)

properties
    surface        % SurfaceMesh carrying the P1 nodes
    order          % Gauss points per triangle: 1, 3, or 7
    projection     % X(nx3) -> X(nx3), snaps a point onto the surface
    nodes6         % quadratic geometry nodes per triangle (nTris x 6 x 3)
    points         % Gauss points on the curved surface ((nTris*order) x 3)
    weights        % curved quadrature weights ((nTris*order) x 1)
    triangleIndex  % source triangle of each point ((nTris*order) x 1)
    basis          % sparse (nPoints x nNodes): P1 basis at each point
end

methods
    function quad = CurvedPanelQuadrature(surface, order, projection)
        arguments
            surface (1,1) SurfaceMesh
            order (1,1) double {mustBeMember(order, [1 3 7])} = 7
            projection (1,1) function_handle = @(X) X
        end
        [bary, wfrac] = triangleGaussRule(order);
        [N, dNdu, dNdv] = CurvedPanelQuadrature.p2Shapes(bary);  % order x 6 each

        tri = surface.tri;
        vtx = surface.vtx;
        nTris = size(tri, 1);
        nNodes = size(vtx, 1);

        quad.surface = surface;
        quad.order = order;
        quad.projection = projection;
        quad.nodes6 = CurvedPanelQuadrature.buildNodes(vtx, tri, projection);

        quad.points = zeros(nTris * order, 3);
        quad.weights = zeros(nTris * order, 1);
        quad.triangleIndex = zeros(nTris * order, 1);

        rows = zeros(nTris * order * 3, 1);
        cols = zeros(nTris * order * 3, 1);
        vals = zeros(nTris * order * 3, 1);
        cursor = 1;
        for t = 1:nTris
            span = (t - 1) * order + (1:order);
            P = squeeze(quad.nodes6(t, :, :));      % 6 x 3 geometry nodes
            Xg = N * P;                             % order x 3 curved points
            Xu = dNdu * P;                          % order x 3 tangents
            Xv = dNdv * P;
            cr = cross(Xu, Xv, 2);                  % order x 3 (2*area normal)
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
        %CENTROIDS Curved panel centroids: the quadratic map at bary [1/3 1/3 1/3].
        N = CurvedPanelQuadrature.p2Shapes([1 1 1] / 3);   % 1 x 6
        nTris = size(quad.nodes6, 1);
        c = zeros(nTris, 3);
        for t = 1:nTris
            c(t, :) = N * squeeze(quad.nodes6(t, :, :));
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
        [~, dNdu, dNdv] = CurvedPanelQuadrature.p2Shapes(bary);
        nTris = size(quad.nodes6, 1);
        Nrm = zeros(quad.nPoints(), 3);
        for t = 1:nTris
            span = (t - 1) * quad.order + (1:quad.order);
            P = squeeze(quad.nodes6(t, :, :));
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

    function nodes6 = buildNodes(vtx, tri, projection)
        %BUILDNODES Quadratic geometry nodes [n1 n2 n3 m12 m23 m31] per triangle.
        %
        % Corners stay on the surface (already meshed there); each edge
        % midpoint is projected onto the surface.  Projection = @(X) X leaves
        % the midpoints straight, so the quadratic map is exactly affine.
        nTris = size(tri, 1);
        nodes6 = zeros(nTris, 6, 3);
        P1 = vtx(tri(:, 1), :);
        P2 = vtx(tri(:, 2), :);
        P3 = vtx(tri(:, 3), :);
        nodes6(:, 1, :) = P1;
        nodes6(:, 2, :) = P2;
        nodes6(:, 3, :) = P3;
        nodes6(:, 4, :) = projection((P1 + P2) / 2);
        nodes6(:, 5, :) = projection((P2 + P3) / 2);
        nodes6(:, 6, :) = projection((P3 + P1) / 2);
    end

    function [N, dNdu, dNdv] = p2Shapes(bary)
        %P2SHAPES Quadratic triangle shapes at barycentric points (u = L2, v = L3).
        %
        % Node order [n1 n2 n3 m12 m23 m31] <-> [L1(2L1-1) L2(2L2-1) L3(2L3-1)
        % 4L1L2 4L2L3 4L3L1].  Returns N, dN/du, dN/dv (each nPts x 6) with
        % L1 = 1-u-v so dL1/du = dL1/dv = -1, dL2/du = 1, dL3/dv = 1.
        L1 = bary(:, 1); L2 = bary(:, 2); L3 = bary(:, 3);
        z = zeros(size(L1));
        N = [L1 .* (2*L1 - 1), L2 .* (2*L2 - 1), L3 .* (2*L3 - 1), ...
             4*L1.*L2, 4*L2.*L3, 4*L3.*L1];
        dNdu = [-(4*L1 - 1), 4*L2 - 1, z, 4*(L1 - L2), 4*L3, -4*L3];
        dNdv = [-(4*L1 - 1), z, 4*L3 - 1, -4*L2, 4*L2, 4*(L1 - L3)];
    end
end
end
