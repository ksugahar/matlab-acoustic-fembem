function sol = femBemCoupledSolve(model, options)
%FEMBEMCOUPLEDSOLVE Johnson-Nedelec coupled FEM/BEM open-boundary solve.
%
% Laplace (default, k = 0):
%
%   -div(c grad u) = f   in the meshed volume (P1 FEM)
%            Delta u = 0   in the unbounded exterior (P1 Galerkin BEM)
%   u and c du/dn continuous on the boundary,  u -> 0 at infinity
%
% Helmholtz ("Wavenumber", k0 > 0) - the ACOUSTIC transmission problem,
% pure MATLAB FEM/BEM:
%
%   (1/rho1) (Delta u + k1^2 u) = -f      in the meshed volume (medium 1)
%    Delta u_e + k0^2 u_e       = 0       in the exterior (medium 0)
%   u = u_e,  (1/rho1) du/dn = (1/rho0) du_e/dn   on the boundary
%   u_e = u_inc + u_s,  u_s outgoing (Sommerfeld via the BEM kernel)
%
% with u_inc = IncidentAmplitude * exp(1i*k0*z), rho0 = 1,
% k1 = InteriorWavenumber (default k0), rho1 = DensityRatio. The coupled
% system is the same Johnson-Nedelec pair with sigma = du_s/dn:
%
%   FEM row :  (1/rho1)(A - k1^2 Mv) u - T' M sigma = F + T' G_inc
%   BIE row :  (1/2 M - K_k) T u + V_k sigma = (1/2 M - K_k) g_inc
%
% where G_inc is the boundary P1 load of du_inc/dn (exact per-triangle
% outward normals, Gauss quadrature) and g_inc the nodal incident trace.
% The exterior scattered representation behind the BIE row is
%
%   u_s(x) = -S_k[sigma](x) + D_k[u_s|_Gamma](x),
%
% exposed through sol.exteriorPotentialAt (the SCATTERED field; add the
% incident wave for totals). For k0 = 0 everything reduces exactly to the
% verified Laplace path (same kind string, same checks).
%
% Analytic gates (tests): unit ball, c = 1, f = 1 (Laplace):
%   u(r) = 1/2 - r^2/6, u_Gamma = 1/3, lambda = -1/3;
% Helmholtz: k1 = k0, rho1 = 1 makes the sphere acoustically INVISIBLE
% (u == u_inc, u_s == 0, exact); k1 ~= k0 is gated by the Anderson fluid-
% sphere series (fluidSphereScattering).

arguments
    model (1,1) FemBemModel
    options.VolumeSource (1,1) double = 1.0
    options.MaterialCoef double = 1
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
    options.Wavenumber (1,1) double {mustBeNonnegative} = 0.0
    options.InteriorWavenumber double = []
    options.DensityRatio (1,1) double {mustBePositive} = 1.0
    options.IncidentAmplitude (1,1) double = 0.0
end

k0 = options.Wavenumber;
k1 = options.InteriorWavenumber;
if isempty(k1)
    k1 = k0;
end
rhor = options.DensityRatio;
amp = options.IncidentAmplitude;
if k0 == 0 && (amp ~= 0 || rhor ~= 1 || k1 ~= 0)
    error("femBemCoupledSolve:laplace", ...
        "IncidentAmplitude / DensityRatio / InteriorWavenumber need Wavenumber > 0.");
end
if k0 > 0 && ~isequal(options.MaterialCoef, 1)
    error("femBemCoupledSolve:material", ...
        "The Helmholtz path uses DensityRatio, not MaterialCoef.");
end

surface = model.surface;
nFem = size(model.mesh.vtx, 1);
nBem = numel(model.mesh.traceNodeIds);

[A, femDetail] = model.h1.stiffness(options.MaterialCoef);
[M, ~] = model.scalarBem.mass();
T = model.trace.matrix;
V = GalerkinSingleLayer(surface, "Wavenumber", k0, ...
    "QuadratureOrder", options.QuadratureOrder);
K = GalerkinDoubleLayer(surface, "Wavenumber", k0, ...
    "QuadratureOrder", options.QuadratureOrder);
quad = V.quadrature;

if k0 > 0
    [Mv, ~] = model.h1.mass(1);
    Afem = (1 / rhor) * (A - k1^2 * Mv);
else
    Afem = A;
end

% P1 load vector of the constant volume source: int f phi_i dV = f*vol/4
F = zeros(nFem, 1);
volumes = femDetail.volumes;
for e = 1:size(model.mesh.tet, 1)
    ids = model.mesh.tet(e, :);
    F(ids) = F(ids) + options.VolumeSource * volumes(e) / 4;
end

% incident plane wave data (zero vectors for amp = 0)
gInc = zeros(nBem, 1);
GInc = zeros(nBem, 1);
if amp ~= 0
    gInc = amp * exp(1i * k0 * surface.vtx(:, 3));
    normals = quad.outwardNormals();
    dnInc = amp * 1i * k0 * normals(:, 3) .* exp(1i * k0 * quad.points(:, 3));
    GInc = quad.weightedBasis().' * dnInc;
end

lhs = [Afem, -T.' * M; (0.5 * M - K.matrix) * T, V.matrix];
rhs = [F + T.' * GInc; (0.5 * M - K.matrix) * gInc];
x = lhs \ rhs;

u = x(1:nFem);
lambda = x(nFem + 1:end);
residual = lhs * x - rhs;
scatteredTrace = T * u - gInc;

sol = struct();
if k0 > 0
    sol.kind = "johnson_nedelec_coupled_fem_bem_helmholtz_solve";
else
    sol.kind = "johnson_nedelec_coupled_fem_bem_solve";
end
sol.policy = "readable_p1_fem_p1_galerkin_bem_open_boundary_teaching_solve";
sol.wavenumber = k0;
sol.interiorWavenumber = k1;
sol.densityRatio = rhor;
sol.incidentAmplitude = amp;
sol.u = u;
sol.lambda = lambda;
sol.trace = T * u;
sol.incidentTrace = gInc;
sol.scatteredTrace = scatteredTrace;
sol.totalExteriorFlux = sum(M * lambda);     % int lambda dS
sol.volumeSourceTotal = options.VolumeSource * sum(volumes);
sol.residualNorm = norm(residual);
sol.quadratureOrder = options.QuadratureOrder;
sol.exteriorPotentialAt = @(points) exteriorPotentialAt( ...
    surface, lambda, scatteredTrace, points, k0, quad);
sol.checks = struct( ...
    "solveResidualSmall", norm(residual) <= 1e-10 * max(1, norm(rhs)), ...
    "fluxBalancesSource", k0 > 0 || ...
        abs(sol.totalExteriorFlux + sol.volumeSourceTotal) ...
        <= 0.05 * max(abs(sol.volumeSourceTotal), eps), ...
    "solutionFieldTypeMatchesKernel", ...
        (k0 == 0) == (isreal(u) && isreal(lambda)));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end


function u = exteriorPotentialAt(surface, lambda, trace, points, k, quad)
%EXTERIORPOTENTIALAT Direct representation u_s = -S_k[lambda] + D_k[trace].
%
% Same split as the operators: analytic Laplace panels for the singular
% parts, smooth low-frequency-stable corrections by quadrature for k > 0
% (so the k -> 0 limit is exactly the verified Laplace evaluation).

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

if k > 0
    parts = HelmholtzKernel(points, quad.points, ...
        "Wavenumber", k, ...
        "SourceWeights", quad.weights, ...
        "SourceNormals", quad.outwardNormals());
    % corrections carry source weights + 1/(4*pi)
    u = u - parts.singleLayerCorrection * (quad.basis * lambda) ...
          + parts.doubleLayerSourceNormalCorrection * (quad.basis * trace);
end
end
