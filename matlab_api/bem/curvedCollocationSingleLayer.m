function op = curvedCollocationSingleLayer(surface, options)
%CURVEDCOLLOCATIONSINGLELAYER P1 nodal-collocation single layer on curved panels.
%
%   op = curvedCollocationSingleLayer(surface, "Wavenumber", k, ...
%           "Projection", CurvedPanelQuadrature.sphereProjection(R));
%   op.matrix          % (nNodes x nNodes) collocation single-layer A_iJ
%   op.potentialMap(x) % (nPts x nNodes) exterior single-layer potential map
%
% Collocation single layer for the exterior Helmholtz/Laplace Dirichlet BVP:
%
%   A_iJ = int_S G_k(x_i, y) phi_J(y) dS(y),   G_k = exp(1i*k*r)/(4*pi*r),
%
% x_i the P1 nodes (mesh vertices, on the surface), phi_J the P1 basis, the
% surface a CURVED (quadratic-isoparametric) panel mesh via CurvedPanelQuadrature.
% Projection = @(X) X gives the flat (straight-panel) reference; a surface
% projection curves the panels -- the ONLY difference between the two, so the
% A/B isolates the O(h^2) faceting error (see curvedVsFlatSoftSphere).
%
% Singular handling is honest and geometry-aware:
%   far panels (x_i not a vertex of the panel)  -> plain curved Gauss
%   incident panels (x_i is a vertex)           -> Duffy-regularised curved
%       quadrature: the reference triangle is collapsed onto the singular
%       vertex, whose Jacobian ~xi cancels the 1/r singularity, integrated
%       on the SAME curved panel geometry so the near field is curved too.
% One-ring near-singular panels use the plain high-order Gauss rule (a
% documented teaching-lane approximation, not a hidden fallback).

arguments
    surface (1,1) SurfaceMesh
    options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 7
    options.Projection (1,1) function_handle = @(X) X
    options.CurveOrder (1,1) double {mustBeMember(options.CurveOrder, [1 2 3])} = 2
    options.GeomNodes double = []          % explicit curved nodes (e.g. from netgen); overrides Projection
    options.DuffyOrder (1,1) double {mustBePositive, mustBeInteger} = 6
end

k = options.Wavenumber;
if isempty(options.GeomNodes)
    quad = CurvedPanelQuadrature(surface, options.QuadratureOrder, ...
        options.Projection, options.CurveOrder);
else
    quad = CurvedPanelQuadrature(surface, options.QuadratureOrder, ...
        "GeomNodes", options.GeomNodes);
end
tri = surface.tri;
colloc = surface.vtx;                 % P1 collocation nodes (on the surface)
nN = size(colloc, 1);
nT = size(tri, 1);
order = options.QuadratureOrder;

% --- far field: every panel by plain curved Gauss (exact for x_i off panel) ---
K = greenKernel(colloc, quad.points, k);            % nN x nGauss
A = (K .* quad.weights.') * quad.basis;             % nN x nN

% --- singular correction: replace each incident panel's plain-Gauss row by a
%     Duffy-regularised curved integral (subtract the wrong part, add the right) -
[gx, gw] = gaussLegendre01(options.DuffyOrder);
for t = 1:nT
    nodes = tri(t, :);
    span = (t - 1) * order + (1:order);
    P = squeeze(quad.geomNodes(t, :, :));           % M x 3 curved geometry nodes
    Kt = quad.weights(span).';                      % 1 x order panel weights
    for p = 1:3
        i = nodes(p);
        xi = colloc(i, :);
        plainRow = (greenKernel(xi, quad.points(span, :), k) .* Kt) ...
            * quad.basis(span, nodes);              % 1 x 3 wrong (singular) part
        duffyRow = duffyPanelRow(P, p, xi, k, gx, gw, quad.curveOrder);   % 1 x 3
        A(i, nodes) = A(i, nodes) - plainRow + duffyRow;
    end
end

op = struct();
if k > 0
    op.kind = "curved_collocation_helmholtz_single_layer";
else
    op.kind = "curved_collocation_laplace_single_layer";
end
op.policy = "p1_nodal_collocation_isoparametric_curved_panel_single_layer";
op.wavenumber = k;
op.matrix = A;
op.quadrature = quad;
op.collocationPoints = colloc;
op.projection = options.Projection;
op.potentialMap = @(points) (greenKernel(points, quad.points, k) ...
    .* quad.weights.') * quad.basis;
end


function G = greenKernel(X, Y, k)
%GREENKERNEL exp(1i*k*r)/(4*pi*r) between rows of X and rows of Y (no toolbox).
D2 = sum(X.^2, 2) + sum(Y.^2, 2).' - 2 * (X * Y.');
D = sqrt(max(D2, 0));
G = exp(1i * k * D) ./ (4 * pi * D);
end


function row = duffyPanelRow(P, pv0, xi, k, gx, gw, curveOrder)
%DUFFYPANELROW int_T G_k(xi,y) [L1 L2 L3](y) dS over a curved panel, singular
% at local vertex pv0.  Duffy: barycentric L(pv0)=1-s, and the opposite edge
% is swept by t in [0,1]; du dv = s ds dt cancels the 1/r vertex singularity.

pv1 = mod(pv0, 3) + 1;
pv2 = mod(pv0 + 1, 3) + 1;
row = zeros(1, 3);
for a = 1:numel(gx)
    s = gx(a);
    for b = 1:numel(gx)
        tt = gx(b);
        L = zeros(1, 3);
        L(pv0) = 1 - s;
        L(pv1) = s * (1 - tt);
        L(pv2) = s * tt;
        [N, dNdu, dNdv] = CurvedPanelQuadrature.lagrangeShapes(curveOrder, L);
        y = N * P;
        cr = cross(dNdu * P, dNdv * P, 2);
        detJ = sqrt(sum(cr.^2, 2));
        r = norm(xi - y);
        G = exp(1i * k * r) / (4 * pi * r);
        w = gw(a) * gw(b) * s * detJ;      % Duffy jac s * surface jac detJ
        row = row + w * G * L;
    end
end
end


function [x, w] = gaussLegendre01(n)
%GAUSSLEGENDRE01 n-point Gauss-Legendre nodes/weights on [0,1] (Golub-Welsch).
beta = (1:n-1) ./ sqrt(4 * (1:n-1).^2 - 1);
Jm = diag(beta, 1) + diag(beta, -1);
[V, Dd] = eig(Jm);
[x, idx] = sort(diag(Dd));
w = (2 * V(1, idx).^2).';
x = 0.5 * (x + 1);                % [-1,1] -> [0,1]
w = 0.5 * w;
end
