classdef SurfaceQuadrature
%SURFACEQUADRATURE Gauss points on the triangles of a SurfaceMesh.
%
%   quad = SurfaceQuadrature(surface, 3);
%   quad.points          % all Gauss points ((nTris*n) x 3)
%   quad.weights         % Gauss weights including triangle areas
%   quad.basis           % sparse (nPoints x nNodes) P1 values at the points
%
% This is the readable counterpart of Gypsilab's dom(mesh, gss): the
% quadrature that Galerkin boundary operators integrate their test (and
% smooth source) factors with.  Supported rules per triangle:
%
%   1  centroid rule            (degree 1)
%   3  interior 3-point rule    (degree 2)
%   7  Dunavant 7-point rule    (degree 5)

properties
    surface        % SurfaceMesh
    order          % points per triangle: 1, 3, or 7
    points         % Gauss points ((nTris*order) x 3)
    weights        % quadrature weights ((nTris*order) x 1), area included
    triangleIndex  % source triangle of each point ((nTris*order) x 1)
    basis          % sparse (nPoints x nNodes): P1 basis at each point
end

methods
    function quad = SurfaceQuadrature(surface, order)
        arguments
            surface (1,1) SurfaceMesh
            order (1,1) double {mustBeMember(order, [1 3 7])} = 3
        end
        [bary, w] = triangleGaussRule(order);
        tri = surface.tri;
        vtx = surface.vtx;
        nTris = size(tri, 1);
        nNodes = size(vtx, 1);
        areas = surface.areas();

        quad.surface = surface;
        quad.order = order;
        quad.points = zeros(nTris * order, 3);
        quad.weights = zeros(nTris * order, 1);
        quad.triangleIndex = zeros(nTris * order, 1);

        rows = zeros(nTris * order * 3, 1);
        cols = zeros(nTris * order * 3, 1);
        vals = zeros(nTris * order * 3, 1);
        cursor = 1;
        for t = 1:nTris
            span = (t - 1) * order + (1:order);
            quad.points(span, :) = bary * vtx(tri(t, :), :);
            quad.weights(span) = w * areas(t);
            quad.triangleIndex(span) = t;
            for g = 1:order
                for k = 1:3
                    rows(cursor) = span(g);
                    cols(cursor) = tri(t, k);
                    vals(cursor) = bary(g, k);   % P1 basis == barycentric
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
end
end


function [bary, w] = triangleGaussRule(order)
%TRIANGLEGAUSSRULE Barycentric points and unit-area weights.

switch order
    case 1
        bary = [1 1 1] / 3;
        w = 1;
    case 3
        bary = [4 1 1; 1 4 1; 1 1 4] / 6;
        w = [1; 1; 1] / 3;
    case 7
        a1 = (6 - sqrt(15)) / 21;
        a2 = (6 + sqrt(15)) / 21;
        w1 = (155 - sqrt(15)) / 1200;
        w2 = (155 + sqrt(15)) / 1200;
        bary = [
            1/3      1/3      1/3
            a1       a1       1 - 2*a1
            a1       1 - 2*a1 a1
            1 - 2*a1 a1       a1
            a2       a2       1 - 2*a2
            a2       1 - 2*a2 a2
            1 - 2*a2 a2       a2];
        w = [9/40; w1; w1; w1; w2; w2; w2];
end
end
