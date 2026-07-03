function sol = rigidScatteringSolve(surface, options)
%RIGIDSCATTERINGSOLVE Sound-hard scattering via the total-field K equation.
%
%   sol = rigidScatteringSolve(surface, "Wavenumber", k);              % direct
%   sol = rigidScatteringSolve(surface, "Wavenumber", k, ...
%       "Method", "chief");                                            % CHIEF
%   sol.trace              % total surface pressure (P1 nodal)
%   u_s = sol.scatteredAt(points);
%
% For a rigid (dp/dn = 0) scatterer the Kirchhoff-Helmholtz representation
% of the TOTAL field loses its single-layer term, so the exterior field is
%
%   u(x) = u_inc(x) + D_k[t](x),        t = total surface trace,
%
% and the boundary limit gives the second-kind equation (in the repo's
% locked outward-normal PV convention; mode-verified on the sphere:
% (1/2 - K_l) t_l = a_l j_l with t_l = a_l (i/x^2)/h_l'(x)):
%
%   (1/2 M - K_k) t = M g_inc.
%
% IRREGULAR FREQUENCIES (taught, then fixed): this equation is singular at
% the interior DIRICHLET eigenvalues of the surface (1/2 - K_l =
% -i k^2 j_l(k) h_l'(k) vanishes where j_l(k R) = 0; unit sphere: first at
% kR = pi). Method "chief" (Schenck 1968) appends interior null-field rows
%
%   D_k[t](x_int) = -u_inc(x_int)
%
% at a few interior points and solves the overdetermined system by least
% squares - the readable classic fix. Burton-Miller (1971) is the
% production alternative (needs the hypersingular operator).

arguments
    surface (1,1) SurfaceMesh
    options.Wavenumber (1,1) double {mustBePositive}
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
    options.Method (1,1) string {mustBeMember(options.Method, ["direct", "chief"])} = "direct"
    options.IncidentAmplitude (1,1) double = 1.0
    options.ChiefPoints double = []
end

k = options.Wavenumber;
amp = options.IncidentAmplitude;
K = GalerkinDoubleLayer(surface, "Wavenumber", k, ...
    "QuadratureOrder", options.QuadratureOrder);
[M, ~] = SurfaceP1Space(surface).mass();
nNodes = size(surface.vtx, 1);

gInc = amp * exp(1i * k * surface.vtx(:, 3));
A = 0.5 * M - K.matrix;
b = M * gInc;

chiefPoints = options.ChiefPoints;
if options.Method == "chief"
    if isempty(chiefPoints)
        % default: centroid plus three jittered interior points (jitter
        % keeps them off the nodal surfaces of the offending interior mode)
        c = mean(surface.vtx, 1);
        span = max(surface.vtx, [], 1) - min(surface.vtx, [], 1);
        chiefPoints = c + 0.15 * span .* ...
            [0 0 0; 0.9 0.3 -0.5; -0.6 0.8 0.4; 0.2 -0.7 0.9];
    end
    C = doubleLayerPotentialMatrix(surface, chiefPoints, k, options.QuadratureOrder);
    cRhs = -amp * exp(1i * k * chiefPoints(:, 3));
    % match the row magnitude of the Galerkin block (M entries ~ area/12)
    w = mean(abs(diag(M)));
    lhs = [A; w * C];
    rhs = [b; w * cRhs];
    t = lhs \ rhs;                      % dense least squares
    residual = lhs * t - rhs;
else
    t = A \ b;
    residual = A * t - b;
end

sol = struct();
sol.kind = "rigid_scattering_total_field_double_layer_solve";
sol.policy = "readable_second_kind_galerkin_p1_bem_teaching_solve";
sol.method = options.Method;
sol.wavenumber = k;
sol.incidentAmplitude = amp;
sol.trace = t;
sol.incidentTrace = gInc;
sol.scatteredTrace = t - gInc;
sol.residualNorm = norm(residual);
sol.conditionNumber = cond(A);
sol.quadratureOrder = options.QuadratureOrder;
sol.chiefPoints = chiefPoints;
sol.operator = K;
sol.surfaceMass = M;
sol.scatteredAt = @(points) ...
    doubleLayerPotentialMatrix(surface, points, k, options.QuadratureOrder) * t;
sol.checks = struct( ...
    "solveResidualSmall", ...
        norm(residual) <= 1e-8 * max(1, norm(b)), ...
    "traceComplex", ~isreal(t));
if all(structfun(@(v) logical(v), sol.checks))
    sol.status = "ok";
else
    sol.status = "needs_attention";
end
end
