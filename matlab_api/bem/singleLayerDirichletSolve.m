function sol = singleLayerDirichletSolve(surface, boundaryValues, options)
%SINGLELAYERDIRICHLETSOLVE Exterior Laplace/Helmholtz Dirichlet BVP via a single layer.
%
%   u(x) = int_S G_k(x,y) q(y) dS(y),   u = g on S,
%   G_k  = exp(1i*k*|x-y|) / (4*pi*|x-y|)     (k = 0: Laplace 1/(4*pi*r))
%
% First-kind boundary integral equation, Galerkin P1 on the boundary:
%
%   V_k q = M g
%
% with V_k the GalerkinSingleLayer matrix (analytic Laplace panels + smooth
% low-frequency-stable Helmholtz correction) and M the boundary P1 surface
% mass. The e^{+ikr} kernel carries the Sommerfeld radiation condition, so
% for k > 0 the same solve is the exterior acoustic Dirichlet problem -
% sound-soft scattering when g = -p_inc on the surface. This is the
% exterior-BEM rung of the cross-validation ladder (stage 4; the Helmholtz
% leg opens the acoustic lane).
%
% For g = 1 on a sphere of radius R the total charge equals the capacitance,
% analytically 4*pi*R at k = 0 (with epsilon0 = 1); for k > 0 sum(M*q) is
% the monopole moment of the density.
%
% Caveat (taught, not hidden): the first-kind V_k equation is singular at
% the interior Dirichlet eigenvalues of the surface (unit sphere: first at
% kR = pi), the classic irregular-frequency limitation addressed later by
% CHIEF / Burton-Miller.
%
%   sol = singleLayerDirichletSolve(surface, ones(nNodes, 1));         % Laplace
%   sol = singleLayerDirichletSolve(surface, g, "Wavenumber", 2.0);    % Helmholtz
%   sol.q              % P1 density (complex for k > 0)
%   sol.totalCharge    % int q dS
%   sol.potentialAt(x) % exterior evaluation, same kernel split as the operator

arguments
    surface (1,1) SurfaceMesh
    boundaryValues double
    options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
end

g = boundaryValues(:);
nNodes = size(surface.vtx, 1);
if numel(g) ~= nNodes
    error("singleLayerDirichletSolve:boundary", ...
        "boundaryValues must have one entry per boundary node (%d).", nNodes);
end

k = options.Wavenumber;
op = GalerkinSingleLayer(surface, "Wavenumber", k, ...
    "QuadratureOrder", options.QuadratureOrder);
[M, ~] = SurfaceP1Space(surface).mass();

rhs = M * g;
q = op.matrix \ rhs;
residualNorm = norm(op.matrix * q - rhs);
totalCharge = sum(M * q);          % int q dS through the exact P1 mass

sol = struct();
if k > 0
    sol.kind = "helmholtz_single_layer_exterior_dirichlet_solve";
else
    sol.kind = "laplace_single_layer_exterior_dirichlet_solve";
end
sol.policy = "readable_first_kind_galerkin_p1_bem_teaching_solve";
sol.wavenumber = k;
sol.q = q;
sol.boundaryValues = g;
sol.totalCharge = totalCharge;
sol.residualNorm = residualNorm;
sol.quadratureOrder = options.QuadratureOrder;
sol.operator = op;
sol.surfaceMass = M;
sol.potentialAt = @(points) potentialAt(surface, q, points, k, op.quadrature);
sol.checks = struct( ...
    "solveResidualSmall", residualNorm <= 1e-10 * max(1, norm(rhs)), ...
    "totalChargeFinite", isfinite(totalCharge), ...
    "densityTypeMatchesKernel", (k == 0) == isreal(q));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end


function u = potentialAt(surface, q, points, k, quad)
%POTENTIALAT Single-layer potential of the solved density at exterior points.
%
% Same split as the operator assembly: the singular Laplace part goes
% through the analytic panel integrals (exact per triangle, robust near the
% surface); for k > 0 the smooth low-frequency-stable Helmholtz correction
% (exp(1i*k*r) - 1)/(4*pi*r) is added by quadrature over the same surface
% Gauss points, so the k -> 0 limit reproduces the Laplace evaluation
% exactly.

u = singleLayerPotentialMatrix(surface, points, k, quad.order) * q;
end
