function sol = curvedSingleLayerDirichletSolve(surface, boundaryValues, options)
%CURVEDSINGLELAYERDIRICHLETSOLVE Exterior Dirichlet BVP on CURVED (isoparametric) panels.
%
%   u(x) = int_S G_k(x,y) q(y) dS(y),   u = g on S,   G_k = exp(1i*k*r)/(4*pi*r)
%
% P1 nodal-collocation single layer (A q = g) on a quadratic-isoparametric
% curved-panel mesh -- the curved-panel counterpart of singleLayerDirichletSolve
% (which is flat Galerkin P1).  The Projection option snaps panel edge nodes
% onto the true surface; Projection = @(X) X is the flat reference, so the same
% call solves both lanes and the difference is purely the O(h^2) faceting error:
% on the sound-soft unit sphere the curved lane is ~10-200x closer to the
% analytic partial-wave series than the flat lane at the SAME mesh (measured,
% testCurvedPanelSphereConvergence).
%
%   proj = CurvedPanelQuadrature.sphereProjection(1.0);
%   g = -exp(1i*k*surface.vtx(:,3));                       % sound-soft plane wave
%   sol = curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k, ...
%             "Projection", proj);
%   sol.q               % P1 density (complex for k > 0)
%   sol.totalCharge     % int q dS  (Laplace g=1 on radius R sphere -> 4*pi*R)
%   sol.potentialAt(x)  % exterior single-layer potential of the density

arguments
    surface (1,1) SurfaceMesh
    boundaryValues double
    options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 7
    options.Projection (1,1) function_handle = @(X) X
    options.CurveOrder (1,1) double {mustBeMember(options.CurveOrder, [1 2 3])} = 2
    options.DuffyOrder (1,1) double {mustBePositive, mustBeInteger} = 6
end

g = boundaryValues(:);
nNodes = size(surface.vtx, 1);
if numel(g) ~= nNodes
    error("curvedSingleLayerDirichletSolve:boundary", ...
        "boundaryValues must have one entry per boundary node (%d).", nNodes);
end

k = options.Wavenumber;
op = curvedCollocationSingleLayer(surface, "Wavenumber", k, ...
    "QuadratureOrder", options.QuadratureOrder, ...
    "Projection", options.Projection, "CurveOrder", options.CurveOrder, ...
    "DuffyOrder", options.DuffyOrder);

q = op.matrix \ g;
residualNorm = norm(op.matrix * q - g);
totalCharge = op.quadrature.weights.' * (op.quadrature.basis * q);

sol = struct();
sol.kind = op.kind;
sol.policy = op.policy;
sol.wavenumber = k;
sol.q = q;
sol.boundaryValues = g;
sol.totalCharge = totalCharge;
sol.residualNorm = residualNorm;
sol.quadratureOrder = options.QuadratureOrder;
sol.operator = op;
sol.projection = options.Projection;
sol.potentialAt = @(points) op.potentialMap(points) * q;
sol.checks = struct( ...
    "solveResidualSmall", residualNorm <= 1e-10 * max(1, norm(g)), ...
    "totalChargeFinite", isfinite(totalCharge), ...
    "densityTypeMatchesKernel", (k == 0) == isreal(q));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end
