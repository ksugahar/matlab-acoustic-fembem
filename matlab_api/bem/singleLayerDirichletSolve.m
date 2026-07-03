function sol = singleLayerDirichletSolve(surface, boundaryValues, options)
%SINGLELAYERDIRICHLETSOLVE Exterior Laplace Dirichlet BVP via a single layer.
%
%   u(x) = int_S G(x,y) q(y) dS(y),   u = g on S,   G = 1/(4*pi*r)
%
% First-kind boundary integral equation, Galerkin P1 on the boundary:
%
%   V q = M g
%
% with V the GalerkinSingleLayer matrix and M the boundary P1 surface mass.
% This is the exterior-BEM rung of the cross-validation ladder (stage 4).
% For g = 1 on a sphere of radius R the total charge equals the
% capacitance, analytically 4*pi*R (with epsilon0 = 1).
%
%   sol = singleLayerDirichletSolve(surface, ones(nNodes, 1));
%   sol.q              % P1 charge density
%   sol.totalCharge    % int q dS - the capacitance for unit voltage
%   sol.potentialAt(x) % readable far-field evaluation

arguments
    surface (1,1) SurfaceMesh
    boundaryValues double
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
end

g = boundaryValues(:);
nNodes = size(surface.vtx, 1);
if numel(g) ~= nNodes
    error("singleLayerDirichletSolve:boundary", ...
        "boundaryValues must have one entry per boundary node (%d).", nNodes);
end

op = GalerkinSingleLayer(surface, "QuadratureOrder", options.QuadratureOrder);
[M, ~] = SurfaceP1Space(surface).mass();

rhs = M * g;
q = op.matrix \ rhs;
residualNorm = norm(op.matrix * q - rhs);
totalCharge = sum(M * q);          % int q dS through the exact P1 mass

sol = struct();
sol.kind = "laplace_single_layer_exterior_dirichlet_solve";
sol.policy = "readable_first_kind_galerkin_p1_bem_teaching_solve";
sol.q = q;
sol.boundaryValues = g;
sol.totalCharge = totalCharge;
sol.residualNorm = residualNorm;
sol.quadratureOrder = options.QuadratureOrder;
sol.operator = op;
sol.surfaceMass = M;
sol.potentialAt = @(points) potentialAt(surface, q, points);
sol.checks = struct( ...
    "solveResidualSmall", residualNorm <= 1e-10 * max(1, norm(rhs)), ...
    "totalChargeFinite", isfinite(totalCharge), ...
    "densityReal", isreal(q));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end


function u = potentialAt(surface, q, points)
%POTENTIALAT Single-layer potential of the solved density at exterior points.
%
% Uses the same analytic panel integrals as the operator, so the evaluation
% stays exact per triangle even close to the surface.

u = zeros(size(points, 1), 1);
tri = surface.tri;
vtx = surface.vtx;
for t = 1:size(tri, 1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), points);
    u = u + I1 * q(tri(t, :));
end
u = u / (4 * pi);
end
