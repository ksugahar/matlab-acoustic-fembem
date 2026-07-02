function sol = laplaceDirichletSolve(model, boundaryValues, options)
%LAPLACEDIRICHLETSOLVE Interior Laplace Dirichlet BVP on P1 tetrahedra.
%
%   -div(c grad u) = 0  in the volume,    u = g  on the whole boundary
%
% This is the first boundary-value-problem rung of the cross-validation
% ladder (VOL_FEM_BEM_COUPLING_DESIGN.md, stage 3): the classic readable
% partition-and-eliminate solve, before any BEM kernel enters.
%
%   K = h1 stiffness;  split dofs into boundary B (the trace nodes) and
%   interior I:
%
%       u_B = g
%       K_II u_I = -K_IB g
%
% For constant c and linear boundary data g = a0 + a.x the P1 solution
% reproduces the linear field EXACTLY (the patch test); the test suite
% locks that identity.
%
%   sol = laplaceDirichletSolve(m, g);
%   sol.u                    % nodal solution on all volume nodes
%   sol.interiorResidualNorm % ||K(I,:) u|| after the solve
%   sol.checks

arguments
    model (1,1) FemBemModel
    boundaryValues double
    options.MaterialCoef double = 1
end

[K, femDetail] = model.h1.stiffness(options.MaterialCoef);
nNodes = size(model.mesh.vtx, 1);
boundaryNodeIds = model.mesh.traceNodeIds(:);
interiorNodeIds = setdiff((1:nNodes).', boundaryNodeIds);

g = boundaryValues(:);
if numel(g) ~= numel(boundaryNodeIds)
    error("laplaceDirichletSolve:boundary", ...
        "boundaryValues must have one entry per boundary trace node (%d).", ...
        numel(boundaryNodeIds));
end

u = zeros(nNodes, 1);
u(boundaryNodeIds) = g;
if isempty(interiorNodeIds)
    solver = "boundary_only_no_interior_unknowns";
else
    Kii = K(interiorNodeIds, interiorNodeIds);
    Kib = K(interiorNodeIds, boundaryNodeIds);
    u(interiorNodeIds) = -Kii \ (Kib * g);
    solver = "backslash_on_spd_interior_block";
end

interiorResidual = K(interiorNodeIds, :) * u;
interiorResidualNorm = norm(interiorResidual);
energy = 0.5 * (u.' * K * u);
boundaryReaction = K(boundaryNodeIds, :) * u;   % discrete Neumann data

sol = struct();
sol.kind = "laplace_dirichlet_interior_solve";
sol.policy = "readable_partition_eliminate_p1_teaching_solve";
sol.u = u;
sol.boundaryNodeIds = boundaryNodeIds;
sol.interiorNodeIds = interiorNodeIds;
sol.boundaryValues = g;
sol.boundaryReaction = boundaryReaction;
sol.interiorResidualNorm = interiorResidualNorm;
sol.energy = energy;
sol.materialCoef = femDetail.materialCoef;
sol.solver = solver;
sol.checks = struct( ...
    "boundaryValuesImposedExactly", isequal(u(boundaryNodeIds), g), ...
    "interiorEquationsSatisfied", interiorResidualNorm <= 1e-10 * max(1, norm(g)), ...
    "energyFiniteNonnegative", isfinite(energy) && energy >= -1e-14, ...
    "reactionBalancesToZero", abs(sum(boundaryReaction)) <= 1e-10 * max(1, norm(boundaryReaction)));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end
