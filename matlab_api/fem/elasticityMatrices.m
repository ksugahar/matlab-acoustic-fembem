function [K, M] = elasticityMatrices(mesh, lambda, mu, solidDensity)
%ELASTICITYMATRICES Vector P1 linear-elasticity stiffness and mass on tets.
%
%   [K, M] = elasticityMatrices(mesh, lambda, mu, rho);
%   Kdyn = K - omega^2 * M;   % time-harmonic dynamic stiffness
%
% Assembles the isotropic linear-elasticity operators for the vector P1
% displacement field on the volume tetrahedra of a VolMesh: the interior
% FEM half of the acoustic fluid-structure interaction (FSI) coupling. The
% displacement is 3 components per node, dof (node n, component c) at
% 3*(n-1)+c, and
%
%   K_ij = int Omega  sigma(u_i) : epsilon(u_j),   sigma = lambda tr(eps) I
%                                                        + 2 mu eps,
%   M_ij = int Omega  rho u_i . u_j.
%
% Per element the strain-displacement matrix B (6x12, Voigt with engineering
% shear) is constant (P1), so K_e = Vol * B' D B with the isotropic
% constitutive D; M_e is the P1 consistent mass kron(mass_tet, I3). Lame
% constants relate to wave speeds by lambda = rho(cL^2 - 2 cT^2),
% mu = rho cT^2.

arguments
    mesh (1,1) VolMesh
    lambda (1,1) double
    mu (1,1) double {mustBeNonnegative}
    solidDensity (1,1) double {mustBePositive}
end

vtx = mesh.vtx;
tet = mesh.tet;
nV = size(vtx, 1);
nT = size(tet, 1);

D = [lambda + 2*mu, lambda, lambda, 0, 0, 0;
     lambda, lambda + 2*mu, lambda, 0, 0, 0;
     lambda, lambda, lambda + 2*mu, 0, 0, 0;
     0, 0, 0, mu, 0, 0;
     0, 0, 0, 0, mu, 0;
     0, 0, 0, 0, 0, mu];
massTet = (ones(4) + eye(4));           % P1 tet consistent-mass pattern (/20)
I3 = eye(3);

rows = zeros(144 * nT, 1);
cols = zeros(144 * nT, 1);
valK = zeros(144 * nT, 1);
valM = zeros(144 * nT, 1);
cursor = 1;
for e = 1:nT
    id = tet(e, :);
    X = vtx(id, :);
    C = [ones(4, 1), X];
    vol = abs(det(C)) / 6;
    if vol <= 0
        error("elasticityMatrices:degenerate", "Degenerate tetrahedron %d.", e);
    end
    Cinv = inv(C);
    grad = Cinv(2:4, :);                  % grad(:, i) = grad phi_i (3 x 1)
    B = zeros(6, 12);
    for i = 1:4
        gx = grad(1, i); gy = grad(2, i); gz = grad(3, i);
        B(:, 3*i-2:3*i) = [gx 0 0; 0 gy 0; 0 0 gz; gy gx 0; 0 gz gy; gz 0 gx];
    end
    Ke = vol * (B.' * D * B);
    Me = solidDensity * vol / 20 * kron(massTet, I3);
    dof = reshape([3*id - 2; 3*id - 1; 3*id], [], 1);
    [I, J] = ndgrid(dof, dof);
    span = cursor:(cursor + 143);
    rows(span) = I(:);
    cols(span) = J(:);
    valK(span) = Ke(:);
    valM(span) = Me(:);
    cursor = cursor + 144;
end

K = sparse(rows, cols, valK, 3*nV, 3*nV);
M = sparse(rows, cols, valM, 3*nV, 3*nV);
end
