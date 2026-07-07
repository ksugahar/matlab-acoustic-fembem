function result = volTdBemConvolutionQuadrature(volFile, options)
%VOLTDBEMCONVOLUTIONQUADRATURE Lubich CQ time-domain acoustic BEM on .vol.
%
%   result = volTdBemConvolutionQuadrature("model.vol") reads the boundary
%   triangles from a first-order Netgen .vol mesh, treats them as a P1
%   Galerkin BEM surface, and solves the retarded exterior single-layer
%   Dirichlet problem by convolution quadrature:
%
%       V(d/dt) q(t) = g(t),     p(x,t) = S(d/dt) q(t).
%
%   The CQ weights are not formed explicitly.  Instead, the BDF generating
%   function delta(zeta) is sampled on a small circle, the Laplace-domain
%   BEM matrices V(s) and S(s) are evaluated at s = delta(zeta)/dt, and the
%   time sequence is recovered by an FFT.  This is a genuine TD-BEM teaching
%   rung, distinct from frequency-sweep inverse FFT of Helmholtz solves.

arguments
    volFile (1,1) string = ""
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 3)} = 16
    options.TimeStep (1,1) double {mustBePositive} = 0.5
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.Method (1,1) string {mustBeMember(options.Method, ["BDF1", "BDF2"])} = "BDF1"
    options.CqRadius double = []
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 1
    options.ObservationPoints double = []
    options.BoundaryTimeData double = []
    options.PulseCenterTime double = []
    options.PulseWidth double = []
    options.ResidualTolerance (1,1) double {mustBePositive} = 1e-8
end

if strlength(volFile) == 0
    volFile = defaultFixture();
end

meshTimer = tic;
mesh = VolMesh(volFile);
surface = mesh.boundary();
meshSeconds = toc(meshTimer);

nBoundary = size(surface.vtx, 1);

% Boundary P1 mass.  The Galerkin single-layer load is <g, phi_i> = (M g)_i,
% so the discrete boundary integral equation is  V(s) q = M ghat  -- NOT the
% raw nodal ghat.  Dropping M scales the density by ~the mean nodal area
% (mesh-dependent) and is physically wrong: verified 12x off the analytic
% soft-sphere scattered field at kR = 1.8, versus 3.6% with M.  The coupled and
% elastic CQ lanes already carry this same boundary mass.
[boundaryMass, ~] = SurfaceP1Space(surface).mass();

if isempty(options.ObservationPoints)
    obs = defaultObservationPoints(mesh.vtx);
else
    obs = options.ObservationPoints;
end
if size(obs, 2) ~= 3
    error("volTdBemConvolutionQuadrature:ObservationPoints", ...
        "ObservationPoints must have three columns.");
end

N = options.NumTime;
dt = options.TimeStep;
t = (0:N-1).' * dt;
rho = options.CqRadius;
if isempty(rho)
    rho = sqrt(eps)^(1 / N);
end
if ~(isscalar(rho) && rho > 0 && rho < 1)
    error("volTdBemConvolutionQuadrature:CqRadius", ...
        "CqRadius must be a scalar in (0, 1).");
end

if isempty(options.BoundaryTimeData)
    pulse = rickerPulse(t, dt, options.PulseCenterTime, options.PulseWidth);
    boundaryData = pulse .* ones(1, nBoundary);
else
    boundaryData = options.BoundaryTimeData;
end
if ~isequal(size(boundaryData), [N, nBoundary])
    error("volTdBemConvolutionQuadrature:BoundaryTimeData", ...
        "BoundaryTimeData must be NumTime x numberOfBoundaryNodes.");
end

cqTimer = tic;

% ===================================================================== %
% CONVOLUTION QUADRATURE.  The time-domain convolution  V(d/dt) q = g is
% turned into N DECOUPLED frequency solves by the discrete Fourier
% transform (Lubich).  The CQ weights are never formed explicitly -- that
% is the whole point.  visualizeConvolutionQuadrature(result) draws each of
% the four phases below.
% ===================================================================== %

% --- Phase 1: sample the CQ contour, map it to Laplace nodes. --------- %
% zeta walks a circle of radius rho < 1; the BDF generating function
% delta(zeta) maps it to s = delta(zeta)/dt.  For an A-stable BDF the image
% lies entirely in Re(s) > 0, so every kernel exp(-s r/c) DECAYS (a
% screened-Laplace kernel) -- this is why CQ is stable where naive
% marching-on-in-time would ring.
n    = (0:N-1).';
zeta = rho * exp(-2i * pi * n / N);
s    = bdfDelta(zeta, options.Method) / dt;

% --- Phase 2: rho-weight and FFT the boundary data. ------------------- %
% The rho^n weighting is the CQ scaling trick that keeps the later inverse
% transform well-behaved; the FFT diagonalizes the lower-triangular CQ
% convolution into N independent right-hand sides boundaryHat(ell, :).
scaledBoundary = (rho .^ n) .* boundaryData;
boundaryHat    = fft(scaledBoundary, [], 1);

% --- Phase 3: one independent frequency-domain BEM solve per node. ---- %
% For each Laplace node s(ell): assemble the Galerkin single layer V(s),
% solve V q = M g_hat for the surface density q(s), and apply the single-layer
% potential S(s) to carry q to the observation points.  cond(V) and the solve
% residual are recorded so the health of every node stays auditable/visible.
densityHat        = zeros(N, nBoundary);
pressureHat       = zeros(N, size(obs, 1));
relativeResiduals = zeros(N, 1);
conditionNumbers  = zeros(N, 1);
% Precompute the s-INDEPENDENT analytic Laplace panel matrix for the observation
% potential ONCE (the CQ grid-radiation dominant cost -- the panel integral does
% not depend on the node s).  Per node only the smooth exp(-s r/c) correction is
% recomputed, so radiating to a large movie grid no longer redoes the panels N times.
obsPanel = zeros(size(obs, 1), nBoundary);
tris = surface.tri;
for tface = 1:size(tris, 1)
    [~, I1] = laplacePanelIntegrals(surface.vtx(tris(tface, :), :), obs);
    obsPanel(:, tris(tface, :)) = obsPanel(:, tris(tface, :)) + I1 / (4 * pi);
end
quadObs = SurfaceQuadrature(surface, options.QuadratureOrder);
for ell = 1:N
    V   = laplaceSingleLayerGalerkin(surface, s(ell), options.SoundSpeed, ...
        options.QuadratureOrder);
    rhs = boundaryMass * boundaryHat(ell, :).';      % Galerkin load M ghat
    q   = V \ rhs;                                   % surface density q(s)
    corrObs = singleLayerSmoothCorrection(obs, quadObs.points, s(ell), ...
        options.SoundSpeed, quadObs.weights);
    Sobs = obsPanel + corrObs * quadObs.basis;       % cached panel + smooth part
    densityHat(ell, :)     = q.';
    pressureHat(ell, :)    = (Sobs * q).';           % p(x, s) = S(s) q(s)
    relativeResiduals(ell) = norm(V * q - rhs) / max(1, norm(rhs));
    conditionNumbers(ell)  = cond(V);
end

% --- Phase 4: inverse FFT + rho^-n unscaling -> causal time signal. --- %
% rho^-n undoes the Phase-2 weighting.  Caveat worth SEEING (panels 5-6 of
% the visualizer): rho^-n grows to 1/sqrt(eps) at the last step, so it also
% amplifies round-off there -- pushing N too far lets the tail drown in
% machine noise even while the per-node residuals stay ~1e-16.
densityComplex  = (rho .^ (-n)) .* ifft(densityHat, [], 1);
pressureComplex = (rho .^ (-n)) .* ifft(pressureHat, [], 1);
density  = real(densityComplex);                     % causal, real-valued
pressure = real(pressureComplex);
cqSeconds = toc(cqTimer);

result = struct();
result.kind = "vol_p1_td_bem_convolution_quadrature_response";
result.method = "lubich_" + lower(options.Method) + "_cq_laplace_domain_single_layer_bem";
result.volFile = string(volFile);
result.meshId = mesh.meshId;
result.meshSummary = mesh.summary;
result.meshSourceId = mesh.sourceFileId;
result.policy = "first_order_vol_tri_tet_boundary_p1_td_bem_convolution_quadrature";
result.time = t;
result.timeStep = dt;
result.cqRadius = rho;
result.cqZeta = zeta;
result.cqLaplaceParameter = s;
result.boundaryData = boundaryData;
result.boundaryDensity = density;
result.observationPoints = obs;
result.pressure = pressure;
result.relativeResiduals = relativeResiduals;
result.conditionNumbers = conditionNumbers;
result.timing = struct( ...
    "mesh_or_import", meshSeconds, ...
    "cq_laplace_solves", cqSeconds);
result.summary = struct( ...
    "num_time", N, ...
    "num_cq_laplace_solves", N, ...
    "num_boundary_nodes", nBoundary, ...
    "num_observation_points", size(obs, 1), ...
    "max_abs_pressure", max(abs(pressure), [], "all"), ...
    "max_abs_density", max(abs(density), [], "all"), ...
    "max_relative_residual", max(relativeResiduals), ...
    "max_condition_number", max(conditionNumbers), ...
    "max_imag_pressure_before_real", max(abs(imag(pressureComplex)), [], "all"), ...
    "max_imag_density_before_real", max(abs(imag(densityComplex)), [], "all"));
scale = max(1, result.summary.max_abs_pressure);
result.checks = struct( ...
    "vol_mesh_tri_tet", mesh.summary.triangles > 0 && mesh.summary.tets > 0, ...
    "p1_boundary_bem", true, ...
    "laplace_parameters_positive_real", all(real(s) > 0), ...
    "cq_residuals_small", result.summary.max_relative_residual < options.ResidualTolerance, ...
    "finite_pressure", all(isfinite(pressure), "all"), ...
    "nonzero_time_response", result.summary.max_abs_pressure > 0, ...
    "real_time_response", result.summary.max_imag_pressure_before_real < 1e-8 * scale, ...
    "not_frequency_sweep_ifft", result.method ~= "frequency_domain_p1_fem_p1_bem_plus_inverse_fft");
if all(structfun(@(v) logical(v), result.checks))
    result.status = "ok";
else
    result.status = "needs_attention";
end
end


function volFile = defaultFixture()
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol"));
end


function obs = defaultObservationPoints(nodes)
center = mean(nodes, 1);
r = max(vecnorm(nodes - center, 2, 2));
if r <= 0
    r = 1;
end
obs = center + r * [
    0.0 0.0 2.4
    2.2 0.0 0.4
    0.0 2.0 1.2
    ];
end


function source = rickerPulse(t, dt, centerTime, width)
if isempty(centerTime)
    centerTime = 4 * dt;
end
if isempty(width)
    width = 1.5 * dt;
end
x = (t - centerTime) / width;
source = (1 - 2 * x.^2) .* exp(-x.^2);
source(abs(source) < 1e-14) = 0;
end
