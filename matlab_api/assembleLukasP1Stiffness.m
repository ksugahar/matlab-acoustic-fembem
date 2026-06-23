function fem = assembleLukasP1Stiffness(model, materialCoef)
%ASSEMBLELUKASP1STIFFNESS Assemble scalar P1 tetra FEM stiffness.
%
% Returns the clean-room counterpart of Lukas H1 grad-grad assembly using the
% same geo.conn_matrix vocabulary. materialCoef is scalar or one value per tet.

arguments
    model (1,1) struct
    materialCoef double = 1
end

nodes = model.lukas.geo.nodes;
tet = model.lukas.geo.conn_matrix;
nNodes = size(nodes, 1);
nTets = size(tet, 1);

if isscalar(materialCoef)
    materialCoef = repmat(materialCoef, nTets, 1);
else
    materialCoef = materialCoef(:);
end
if numel(materialCoef) ~= nTets
    error("assembleLukasP1Stiffness:material", "materialCoef must be scalar or one value per tetrahedron.");
end

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

fem = struct();
fem.family = "P1";
fem.cell = "tetrahedron";
fem.stiffness = sparse(ii, jj, vv, nNodes, nNodes);
fem.localStiffness = local;
fem.gradients = gradients;
fem.volumes = volumes;
fem.materialCoef = materialCoef;
end


function [Ke, gradPhi, volume] = localP1TetStiffness(X, materialCoef)
%LOCALP1TETSTIFFNESS Exact P1 stiffness for one tetrahedron.

D = [ones(4, 1), X];
detD = det(D);
volume = abs(detD) / 6;
if volume <= eps
    error("assembleLukasP1Stiffness:degenerate", "Degenerate tetrahedron with near-zero volume.");
end

coeff = D \ eye(4);
gradPhi = coeff(2:4, :).';
Ke = materialCoef * volume * (gradPhi * gradPhi.');
end
