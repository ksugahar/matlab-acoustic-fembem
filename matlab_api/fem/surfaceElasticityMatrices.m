function [K, M] = surfaceElasticityMatrices(surface, youngs, poisson, thickness, density)
%SURFACEELASTICITYMATRICES P1 membrane (plane-stress CST) elasticity on a surface.
%
%   [K, M] = surfaceElasticityMatrices(surface, E, nu, thickness, rho);
%
% The SURFACE-ELEMENT counterpart of elasticityMatrices (which is volume tets):
% a first-order plane-stress constant-strain-triangle (CST) membrane on the
% boundary triangles, embedded in 3D.  Displacement is 3 components per surface
% node, dof (node n, component c) at 3*(n-1)+c -- the SAME layout as
% elasticityMatrices and the interface coupling G, so a shell-FEM / acoustic-BEM
% FSI reuses the identical coupling block with this K,M in place of the volume
% ones.
%
%   K_ij = int_S  t sigma(u_i) : epsilon(u_j),   plane stress,
%   M_ij = int_S  rho t u_i . u_j.
%
% Each element resists only IN-PLANE (membrane) strain; on a CLOSED CURVED
% surface the assembled stiffness has exactly the 6 rigid-body zero modes
% (a node's normal motion stretches its neighbours), so the rigid-body patch
% test ||K u_rigid|| ~ 0 holds -- the self-consistency this operator is here to
% demonstrate for the surface-element variant.
%
% SCOPE: this is the membrane demonstrator for the surface-element FSI variant's
% self-consistency.  A PRODUCTION shell scatterer additionally needs BENDING
% (Kirchhoff / Reissner-Mindlin) for its normal-pressure response, and -- when
% the fluid wets BOTH faces -- an open-surface BEM carrying the pressure jump
% [p] = p+ - p- through the hypersingular operator W (the full Calderon
% projector [[1/2-K, V],[W, 1/2+K']]).  The membrane K here is enough to
% exercise the coupling harness, not to replace a shell solver.

arguments
    surface (1,1) SurfaceMesh
    youngs (1,1) double {mustBePositive} = 1.0
    poisson (1,1) double {mustBeGreaterThan(poisson,-1), mustBeLessThan(poisson,0.5)} = 0.3
    thickness (1,1) double {mustBePositive} = 0.05
    density (1,1) double {mustBePositive} = 1.0
end

vtx = surface.vtx;
tri = surface.tri;
nV = size(vtx, 1);
nT = size(tri, 1);

Dps = youngs / (1 - poisson^2) * ...
    [1, poisson, 0; poisson, 1, 0; 0, 0, (1 - poisson) / 2];   % plane stress

rows = zeros(81 * nT, 1);
cols = zeros(81 * nT, 1);
valK = zeros(81 * nT, 1);
valM = zeros(81 * nT, 1);
cursor = 1;
for e = 1:nT
    id = tri(e, :);
    X = vtx(id, :);
    nrm = cross(X(2, :) - X(1, :), X(3, :) - X(1, :));
    area = 0.5 * norm(nrm);
    if area <= eps
        error("surfaceElasticityMatrices:degenerate", "Degenerate triangle %d.", e);
    end
    nrm = nrm / norm(nrm);
    e1 = (X(2, :) - X(1, :)); e1 = e1 / norm(e1);
    e2 = cross(nrm, e1);                          % in-plane orthonormal frame

    P = [(X - X(1, :)) * e1.', (X - X(1, :)) * e2.'];   % local 2D coords (3x2)
    xl = P(:, 1); yl = P(:, 2);
    b = [yl(2) - yl(3); yl(3) - yl(1); yl(1) - yl(2)];
    c = [xl(3) - xl(2); xl(1) - xl(3); xl(2) - xl(1)];
    B = zeros(3, 6);
    for i = 1:3
        B(:, 2*i-1:2*i) = [b(i), 0; 0, c(i); c(i), b(i)] / (2 * area);
    end
    Kloc = thickness * area * (B.' * Dps * B);   % 6x6 in local (e1,e2) node dofs

    T = zeros(6, 9);                             % local 2D <- global 3D nodal dofs
    for i = 1:3
        T(2*i-1, 3*i-2:3*i) = e1;
        T(2*i,   3*i-2:3*i) = e2;
    end
    Ke = T.' * Kloc * T;                         % 9x9 in 3D nodal dofs

    massTri = area / 12 * [2 1 1; 1 2 1; 1 1 2];
    Me = density * thickness * kron(massTri, eye(3));

    dof = reshape([3*id - 2; 3*id - 1; 3*id], [], 1);
    [I, J] = ndgrid(dof, dof);
    span = cursor:(cursor + 80);
    rows(span) = I(:);
    cols(span) = J(:);
    valK(span) = Ke(:);
    valM(span) = Me(:);
    cursor = cursor + 81;
end

K = sparse(rows, cols, valK, 3*nV, 3*nV);
M = sparse(rows, cols, valM, 3*nV, 3*nV);
end
