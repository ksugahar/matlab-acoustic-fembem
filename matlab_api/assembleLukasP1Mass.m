function fem = assembleLukasP1Mass(model, materialCoef)
%ASSEMBLELUKASP1MASS Assemble scalar P1 tetra FEM mass matrix.

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
    error("assembleLukasP1Mass:material", "materialCoef must be scalar or one value per tetrahedron.");
end

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

fem = struct();
fem.family = "P1";
fem.cell = "tetrahedron";
fem.mass = sparse(ii, jj, vv, nNodes, nNodes);
fem.localMass = local;
fem.volumes = volumes;
fem.materialCoef = materialCoef;
end
