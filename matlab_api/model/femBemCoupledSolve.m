function sol = femBemCoupledSolve(model, options)
%FEMBEMCOUPLEDSOLVE Johnson-Nedelec coupled FEM/BEM open-boundary solve.
%
%   -div(c grad u) = f   in the meshed volume (P1 FEM)
%             Delta u = 0   in the unbounded exterior (P1 Galerkin BEM)
%   u and c du/dn continuous on the boundary,  u -> 0 at infinity
%
% This is the final rung of the cross-validation ladder (stage 5). With
% lambda = du_e/dn (outward normal flux of the exterior field, equal to
% c du/dn of the interior field by transmission), the coupled system is
% the classic Johnson-Nedelec pair:
%
%   FEM row :  A u - T' M lambda = F
%   BIE row :  (1/2 M - K) T u + V lambda = 0
%
% where T is the one-hot trace, M the boundary P1 mass, V the Galerkin
% single layer, and K the Galerkin double layer (outward normal, principal
% value; the +1/2 jump is explicit in the BIE row). The exterior direct
% representation behind the BIE row is
%
%   u_e(x) = -S[lambda](x) + D[u_Gamma](x),
%
% exposed through sol.exteriorPotentialAt for far-field checks.
%
% Analytic gate (tests): unit ball, c = 1, f = 1 gives
%   u(r) = 1/2 - r^2/6,   u_Gamma = 1/3,   lambda = -1/3.

arguments
    model (1,1) FemBemModel
    options.VolumeSource (1,1) double = 1.0
    options.MaterialCoef double = 1
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
end

surface = model.surface;
nFem = size(model.mesh.vtx, 1);
nBem = numel(model.mesh.traceNodeIds);

[A, femDetail] = model.h1.stiffness(options.MaterialCoef);
[M, ~] = model.scalarBem.mass();
T = model.trace.matrix;
V = GalerkinSingleLayer(surface, "QuadratureOrder", options.QuadratureOrder);
K = GalerkinDoubleLayer(surface, "QuadratureOrder", options.QuadratureOrder);

% P1 load vector of the constant volume source: int f phi_i dV = f*vol/4
F = zeros(nFem, 1);
volumes = femDetail.volumes;
for e = 1:size(model.mesh.tet, 1)
    ids = model.mesh.tet(e, :);
    F(ids) = F(ids) + options.VolumeSource * volumes(e) / 4;
end

lhs = [A, -T.' * M; (0.5 * M - K.matrix) * T, V.matrix];
rhs = [F; zeros(nBem, 1)];
x = lhs \ rhs;

u = x(1:nFem);
lambda = x(nFem + 1:end);
residual = lhs * x - rhs;

sol = struct();
sol.kind = "johnson_nedelec_coupled_fem_bem_solve";
sol.policy = "readable_p1_fem_p1_galerkin_bem_open_boundary_teaching_solve";
sol.u = u;
sol.lambda = lambda;
sol.trace = T * u;
sol.totalExteriorFlux = sum(M * lambda);     % int lambda dS
sol.volumeSourceTotal = options.VolumeSource * sum(volumes);
sol.residualNorm = norm(residual);
sol.quadratureOrder = options.QuadratureOrder;
sol.exteriorPotentialAt = @(points) exteriorPotentialAt(surface, lambda, T * u, points);
sol.checks = struct( ...
    "solveResidualSmall", norm(residual) <= 1e-10 * max(1, norm(rhs)), ...
    "fluxBalancesSource", abs(sol.totalExteriorFlux + sol.volumeSourceTotal) ...
        <= 0.05 * max(abs(sol.volumeSourceTotal), eps), ...
    "solutionReal", isreal(u) && isreal(lambda));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end


function u = exteriorPotentialAt(surface, lambda, trace, points)
%EXTERIORPOTENTIALAT Direct representation u_e = -S[lambda] + D[trace].

signs = surface.orientation.triangleOrientationSignsToOutward(:);
u = zeros(size(points, 1), 1);
tri = surface.tri;
vtx = surface.vtx;
for t = 1:size(tri, 1)
    Vt = vtx(tri(t, :), :);
    [~, I1] = laplacePanelIntegrals(Vt, points);
    [~, J1] = laplaceDoubleLayerPanelIntegrals(Vt, points);
    u = u - I1 * lambda(tri(t, :)) / (4 * pi) ...
          + signs(t) * (J1 * trace(tri(t, :))) / (4 * pi);
end
end
