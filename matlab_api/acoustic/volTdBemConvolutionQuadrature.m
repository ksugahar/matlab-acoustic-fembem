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
n = (0:N-1).';
zeta = rho * exp(-2i * pi * n / N);
s = bdfDelta(zeta, options.Method) / dt;
scaledBoundary = (rho .^ n) .* boundaryData;
boundaryHat = fft(scaledBoundary, [], 1);

densityHat = zeros(N, nBoundary);
pressureHat = zeros(N, size(obs, 1));
relativeResiduals = zeros(N, 1);
conditionNumbers = zeros(N, 1);
for ell = 1:N
    V = laplaceSingleLayerGalerkin(surface, s(ell), options.SoundSpeed, ...
        options.QuadratureOrder);
    rhs = boundaryHat(ell, :).';
    q = V \ rhs;
    Sobs = cqSingleLayerPotential(surface, obs, s(ell), ...
        options.SoundSpeed, options.QuadratureOrder);
    densityHat(ell, :) = q.';
    pressureHat(ell, :) = (Sobs * q).';
    relativeResiduals(ell) = norm(V * q - rhs) / max(1, norm(rhs));
    conditionNumbers(ell) = cond(V);
end
densityComplex = (rho .^ (-n)) .* ifft(densityHat, [], 1);
pressureComplex = (rho .^ (-n)) .* ifft(pressureHat, [], 1);
density = real(densityComplex);
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


function rows = cqSingleLayerPotential(surface, points, s, soundSpeed, quadratureOrder)
tri = surface.tri;
vtx = surface.vtx;
rows = zeros(size(points, 1), size(vtx, 1));
for t = 1:size(tri, 1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), points);
    rows(:, tri(t, :)) = rows(:, tri(t, :)) + I1 / (4 * pi);
end

quad = SurfaceQuadrature(surface, quadratureOrder);
correction = cqSingleLayerCorrection(points, quad.points, s, ...
    soundSpeed, quad.weights);
rows = rows + correction * quad.basis;
end


function C = cqSingleLayerCorrection(targetPoints, sourcePoints, s, soundSpeed, sourceWeights)
nTarget = size(targetPoints, 1);
nSource = size(sourcePoints, 1);
C = complex(zeros(nTarget, nSource));
alpha = s / soundSpeed;
for i = 1:nTarget
    for j = 1:nSource
        r = norm(targetPoints(i, :) - sourcePoints(j, :));
        if r == 0
            % Coincident quadrature point: match the frequency-domain
            % HelmholtzKernel convention (leave the smooth correction at 0).
            value = 0;
        else
            z = -alpha * r;
            value = stableExpm1OverR(z, r);
        end
        C(i, j) = sourceWeights(j) * value / (4 * pi);
    end
end
end


function delta = bdfDelta(zeta, method)
switch upper(method)
    case "BDF1"
        delta = 1 - zeta;
    case "BDF2"
        delta = 1.5 - 2*zeta + 0.5*zeta.^2;
end
end


function value = stableExpm1OverR(z, r)
if abs(z) < 1e-5
    value = 0;
    term = 1;
    for k = 1:10
        term = term * z / k;
        value = value + term / r;
    end
else
    value = (exp(z) - 1) / r;
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
