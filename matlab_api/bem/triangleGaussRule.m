function [bary, w] = triangleGaussRule(order)
%TRIANGLEGAUSSRULE Barycentric Gauss points and unit-sum weights on a triangle.
%
%   [bary, w] = triangleGaussRule(order)
%   bary : order x 3 barycentric coordinates (L1, L2, L3) of the Gauss points
%   w    : order x 1 weights summing to 1 (fractions of the triangle)
%
% The physical weight of a point is w * (area of the triangle) for a flat
% panel, or w * 0.5 * |x_u x x_v| for a curved (isoparametric) panel.  Shared
% by SurfaceQuadrature (flat) and CurvedPanelQuadrature (curved) so the two
% integrate the SAME reference rule and differ only in geometry.  Rules:
%
%   1  centroid rule            (degree 1)
%   3  interior 3-point rule    (degree 2)
%   7  Dunavant 7-point rule    (degree 5)

arguments
    order (1,1) double {mustBeMember(order, [1 3 7])} = 3
end

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
