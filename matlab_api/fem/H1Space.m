classdef H1Space
%H1SPACE First-order H1 (P1) space on the tetrahedra of a VolMesh.
%
%   space = H1Space(mesh);
%   [K, detail] = space.stiffness(materialCoef);   % int c grad(u).grad(v) dx
%   [M, detail] = space.mass(materialCoef);        % int c u v dx
%
% The dofs are the volume nodes, so ndof equals size(mesh.vtx, 1). The
% assembly is the classic per-tetrahedron loop kept explicit for teaching:
% barycentric gradients from the 4x4 node matrix, exact P1 integrals, and a
% sparse triplet accumulation.

properties (Constant)
    family = "Lagrange"
    order = 1
    cell = "tetrahedron"
    basis = "P1"
end

properties
    mesh   % VolMesh carrying vtx/tet
end

methods
    function space = H1Space(mesh)
        arguments
            mesh (1,1) VolMesh
        end
        space.mesh = mesh;
    end

    function n = ndof(space)
        n = size(space.mesh.vtx, 1);
    end

    function ids = traceNodeIds(space)
        %TRACENODEIDS Volume node ids that carry the boundary trace.
        ids = space.mesh.traceNodeIds;
    end

    function [K, detail] = stiffness(space, materialCoef)
        %STIFFNESS Assemble int c grad(u).grad(v) dx on P1 tetrahedra.
        arguments
            space (1,1) H1Space
            materialCoef double = 1
        end
        nodes = space.mesh.vtx;
        tet = space.mesh.tet;
        nNodes = size(nodes, 1);
        nTets = size(tet, 1);
        materialCoef = perTetCoefficient(materialCoef, nTets);

        ii = zeros(16 * nTets, 1);
        jj = zeros(16 * nTets, 1);
        vv = zeros(16 * nTets, 1);
        local = zeros(4, 4, nTets);
        gradients = zeros(4, 3, nTets);
        volumes = zeros(nTets, 1);

        cursor = 1;
        for e = 1:nTets
            ids = tet(e, :);
            X = nodes(ids, :);
            [Ke, gradPhi, volume] = localP1TetStiffness(X, materialCoef(e));
            local(:, :, e) = Ke;
            gradients(:, :, e) = gradPhi;
            volumes(e) = volume;

            [I, J] = ndgrid(ids, ids);
            span = cursor:(cursor + 15);
            ii(span) = I(:);
            jj(span) = J(:);
            vv(span) = Ke(:);
            cursor = cursor + 16;
        end

        K = sparse(ii, jj, vv, nNodes, nNodes);
        detail = struct( ...
            "family", space.basis, ...
            "cell", space.cell, ...
            "stiffness", K, ...
            "localStiffness", local, ...
            "gradients", gradients, ...
            "volumes", volumes, ...
            "materialCoef", materialCoef);
    end

    function [M, detail] = mass(space, materialCoef)
        %MASS Assemble int c u v dx with the exact P1 tetra mass rule.
        arguments
            space (1,1) H1Space
            materialCoef double = 1
        end
        nodes = space.mesh.vtx;
        tet = space.mesh.tet;
        nNodes = size(nodes, 1);
        nTets = size(tet, 1);
        materialCoef = perTetCoefficient(materialCoef, nTets);

        ii = zeros(16 * nTets, 1);
        jj = zeros(16 * nTets, 1);
        vv = zeros(16 * nTets, 1);
        local = zeros(4, 4, nTets);
        volumes = zeros(nTets, 1);

        cursor = 1;
        for e = 1:nTets
            ids = tet(e, :);
            X = nodes(ids, :);
            volume = abs(det([ones(4, 1), X])) / 6;
            Me = materialCoef(e) * volume / 20 * (ones(4, 4) + eye(4));
            local(:, :, e) = Me;
            volumes(e) = volume;

            [I, J] = ndgrid(ids, ids);
            span = cursor:(cursor + 15);
            ii(span) = I(:);
            jj(span) = J(:);
            vv(span) = Me(:);
            cursor = cursor + 16;
        end

        M = sparse(ii, jj, vv, nNodes, nNodes);
        detail = struct( ...
            "family", space.basis, ...
            "cell", space.cell, ...
            "mass", M, ...
            "localMass", local, ...
            "volumes", volumes, ...
            "materialCoef", materialCoef);
    end
end
end


function coef = perTetCoefficient(materialCoef, nTets)
%PERTETCOEFFICIENT Expand a scalar coefficient to one value per tetrahedron.

if isscalar(materialCoef)
    coef = repmat(materialCoef, nTets, 1);
else
    coef = materialCoef(:);
end
if numel(coef) ~= nTets
    error("H1Space:material", ...
        "materialCoef must be scalar or one value per tetrahedron.");
end
end


function [Ke, gradPhi, volume] = localP1TetStiffness(X, materialCoef)
%LOCALP1TETSTIFFNESS Exact P1 stiffness for one tetrahedron.

D = [ones(4, 1), X];
detD = det(D);
volume = abs(detD) / 6;
if volume <= eps
    error("H1Space:degenerate", "Degenerate tetrahedron with near-zero volume.");
end

coeff = D \ eye(4);
gradPhi = coeff(2:4, :).';
Ke = materialCoef * volume * (gradPhi * gradPhi.');
end
