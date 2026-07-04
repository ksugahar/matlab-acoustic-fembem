function result = volFemBemCoupledConvolutionQuadrature(volFile, options)
%VOLFEMBEMCOUPLEDCONVOLUTIONQUADRATURE .vol volume-FEM / BEM coupled CQ.
%
%   result = volFemBemCoupledConvolutionQuadrature("model.vol") reads a
%   first-order Netgen .vol mesh, builds H1/P1 tetrahedral volume FEM and
%   P1 boundary BEM spaces, and solves the Johnson-Nedelec / Calderon
%   coupled Laplace-domain system at Lubich CQ points:
%
%       [ A + (s/c1)^2 Mv      -T' Mb ] [ uhat ] = [ Fhat(t) ]
%       [ (1/2 Mb - K(s)) T     V(s)  ] [ qhat ]   [   0     ]
%
%   where V(s) and K(s) are retarded single- and double-layer Galerkin
%   matrices with kernel exp(-s r/c0)/(4*pi*r).  The exterior field is
%   evaluated by the matching representation -S(s)q + D(s)Tu.  The inverse
%   CQ FFT recovers interior pressure u(t), boundary flux density q(t), and
%   exterior pressure p(x,t).  CouplingForm="SingleLayerTeaching" keeps the
%   earlier one-equation rung for regression and comparison.

arguments
    volFile (1,1) string = ""
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 3)} = 16
    options.TimeStep (1,1) double {mustBePositive} = 0.5
    options.ExteriorSoundSpeed (1,1) double {mustBePositive} = 1.0
    options.InteriorSoundSpeed (1,1) double {mustBePositive} = 1.0
    options.Method (1,1) string {mustBeMember(options.Method, ["BDF1", "BDF2"])} = "BDF1"
    options.CouplingForm (1,1) string {mustBeMember(options.CouplingForm, ...
        ["JohnsonNedelec", "SingleLayerTeaching"])} = "JohnsonNedelec"
    options.CqRadius double = []
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 1
    options.ObservationPoints double = []
    options.VolumeSourceTimeData double = []
    options.PulseCenterTime double = []
    options.PulseWidth double = []
    options.ResidualTolerance (1,1) double {mustBePositive} = 1e-8
end

if strlength(volFile) == 0
    volFile = defaultFixture();
end

meshTimer = tic;
model = FemBemModel(volFile);
surface = model.surface;
nVolume = size(model.mesh.vtx, 1);
nBoundary = size(surface.vtx, 1);
[A, ~] = model.h1.stiffness(1);
[Mv, ~] = model.h1.mass(1);
[Mb, ~] = SurfaceP1Space(surface).mass();
T = model.trace.matrix;
meshSeconds = toc(meshTimer);

if isempty(options.ObservationPoints)
    obs = defaultObservationPoints(model.mesh.vtx);
else
    obs = options.ObservationPoints;
end
if size(obs, 2) ~= 3
    error("volFemBemCoupledConvolutionQuadrature:ObservationPoints", ...
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
    error("volFemBemCoupledConvolutionQuadrature:CqRadius", ...
        "CqRadius must be a scalar in (0, 1).");
end

if isempty(options.VolumeSourceTimeData)
    source = rickerPulse(t, dt, options.PulseCenterTime, options.PulseWidth);
else
    source = options.VolumeSourceTimeData(:);
end
if numel(source) ~= N
    error("volFemBemCoupledConvolutionQuadrature:VolumeSourceTimeData", ...
        "VolumeSourceTimeData must have NumTime entries.");
end

cqTimer = tic;
n = (0:N-1).';
zeta = rho * exp(-2i * pi * n / N);
s = bdfDelta(zeta, options.Method) / dt;
sourceHat = fft((rho .^ n) .* source);
unitLoad = Mv * ones(nVolume, 1);
useJohnsonNedelec = options.CouplingForm == "JohnsonNedelec";

Uhat = zeros(N, nVolume);
Qhat = zeros(N, nBoundary);
Phat = zeros(N, size(obs, 1));
relativeResiduals = zeros(N, 1);
conditionNumbers = zeros(N, 1);
doubleLayerNorms = zeros(N, 1);
for ell = 1:N
    alpha = s(ell) / options.InteriorSoundSpeed;
    Kdyn = A + alpha^2 * Mv;
    V = cqSingleLayerGalerkin(surface, s(ell), options.ExteriorSoundSpeed, ...
        options.QuadratureOrder);
    Sobs = cqSingleLayerPotential(surface, obs, s(ell), ...
        options.ExteriorSoundSpeed, options.QuadratureOrder);
    if useJohnsonNedelec
        K = cqDoubleLayerGalerkin(surface, s(ell), options.ExteriorSoundSpeed, ...
            options.QuadratureOrder);
        Dobs = cqDoubleLayerPotential(surface, obs, s(ell), ...
            options.ExteriorSoundSpeed, options.QuadratureOrder);
        couplingRow = (0.5 * Mb - K) * T;
        bemBlock = V;
        fieldRow = @(u, q) -Sobs * q + Dobs * (T * u);
        doubleLayerNorms(ell) = norm(K, "fro");
    else
        couplingRow = Mb * T;
        bemBlock = -V;
        fieldRow = @(u, q) Sobs * q;
        doubleLayerNorms(ell) = 0.0;
    end
    lhs = [Kdyn, -T.' * Mb; couplingRow, bemBlock];
    rhs = [unitLoad * sourceHat(ell); zeros(nBoundary, 1)];
    x = lhs \ rhs;
    u = x(1:nVolume);
    q = x(nVolume + (1:nBoundary));
    Uhat(ell, :) = u.';
    Qhat(ell, :) = q.';
    Phat(ell, :) = fieldRow(u, q).';
    relativeResiduals(ell) = norm(lhs * x - rhs) / max(1, norm(rhs));
    conditionNumbers(ell) = cond(full(lhs));
end
interiorComplex = (rho .^ (-n)) .* ifft(Uhat, [], 1);
densityComplex = (rho .^ (-n)) .* ifft(Qhat, [], 1);
pressureComplex = (rho .^ (-n)) .* ifft(Phat, [], 1);
interior = real(interiorComplex);
density = real(densityComplex);
pressure = real(pressureComplex);
cqSeconds = toc(cqTimer);

result = struct();
result.kind = "vol_p1_volume_fem_boundary_bem_coupled_cq_response";
if useJohnsonNedelec
    methodTail = "_johnson_nedelec_calderon_bem_coupling";
else
    methodTail = "_single_layer_teaching_bem_coupling";
end
result.method = "lubich_" + lower(options.Method) + "_cq_volume_p1_fem" + methodTail;
result.couplingForm = options.CouplingForm;
result.volFile = string(volFile);
result.meshId = model.mesh.meshId;
result.meshSummary = model.mesh.summary;
result.meshSourceId = model.mesh.sourceFileId;
result.policy = "first_order_vol_tri_tet_volume_p1_fem_boundary_p1_bem_cq_coupling";
result.time = t;
result.timeStep = dt;
result.cqRadius = rho;
result.cqZeta = zeta;
result.cqLaplaceParameter = s;
result.volumeSourceTime = source;
result.interiorPressure = interior;
result.boundaryDensity = density;
result.observationPoints = obs;
result.exteriorPressure = pressure;
result.relativeResiduals = relativeResiduals;
result.conditionNumbers = conditionNumbers;
result.doubleLayerNorms = doubleLayerNorms;
result.timing = struct( ...
    "mesh_or_import", meshSeconds, ...
    "coupled_cq_laplace_solves", cqSeconds);

scaleP = max(abs(pressure), [], "all");
scaleU = max(abs(interior), [], "all");
result.summary = struct( ...
    "num_time", N, ...
    "num_coupled_cq_laplace_solves", N, ...
    "num_volume_nodes", nVolume, ...
    "num_boundary_nodes", nBoundary, ...
    "num_interior_only_nodes", nVolume - nBoundary, ...
    "num_observation_points", size(obs, 1), ...
    "max_abs_interior_pressure", scaleU, ...
    "max_abs_boundary_density", max(abs(density), [], "all"), ...
    "max_abs_exterior_pressure", scaleP, ...
    "max_relative_residual", max(relativeResiduals), ...
    "max_condition_number", max(conditionNumbers), ...
    "max_double_layer_frobenius_norm", max(doubleLayerNorms), ...
    "max_imag_interior_before_real", max(abs(imag(interiorComplex)), [], "all"), ...
    "max_imag_exterior_before_real", max(abs(imag(pressureComplex)), [], "all"));
result.checks = struct( ...
    "vol_mesh_tri_tet", model.mesh.summary.triangles > 0 && model.mesh.summary.tets > 0, ...
    "p1_volume_fem", model.h1.order == 1 && model.h1.cell == "tetrahedron", ...
    "p1_boundary_bem", model.scalarBem.order == 1 && model.scalarBem.cell == "triangle", ...
    "laplace_parameters_positive_real", all(real(s) > 0), ...
    "coupled_residuals_small", result.summary.max_relative_residual < options.ResidualTolerance, ...
    "finite_interior", all(isfinite(interior), "all"), ...
    "finite_exterior", all(isfinite(pressure), "all"), ...
    "nonzero_interior_response", scaleU > 0, ...
    "nonzero_exterior_response", scaleP > 0, ...
    "johnson_nedelec_calderon_form", useJohnsonNedelec, ...
    "double_layer_k_included", useJohnsonNedelec && max(doubleLayerNorms) > 0, ...
    "single_layer_teaching_form_available", true, ...
    "real_interior_response", result.summary.max_imag_interior_before_real < 1e-8 * max(1, scaleU), ...
    "real_exterior_response", result.summary.max_imag_exterior_before_real < 1e-8 * max(1, scaleP), ...
    "has_volume_fem_unknowns", nVolume > 0, ...
    "not_exterior_only_td_bem", true);
if all(structfun(@(v) logical(v), result.checks))
    result.status = "ok";
else
    result.status = "needs_attention";
end
end


function K = cqDoubleLayerGalerkin(surface, s, soundSpeed, quadratureOrder)
quad = SurfaceQuadrature(surface, quadratureOrder);
signs = surface.orientation.triangleOrientationSignsToOutward(:);
if any(signs == 0)
    error("volFemBemCoupledConvolutionQuadrature:orientation", ...
        "Surface orientation is unknown for %d triangle(s); cannot assemble K(s).", ...
        sum(signs == 0));
end

nGauss = quad.nPoints();
nNodes = size(surface.vtx, 1);
tri = surface.tri;
vtx = surface.vtx;

P = complex(zeros(nGauss, nNodes));
for t = 1:size(tri, 1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t, :), :), quad.points);
    P(:, tri(t, :)) = P(:, tri(t, :)) + signs(t) * J1;
end
Bw = quad.weightedBasis();
K = Bw.' * P / (4 * pi);

correction = cqDoubleLayerCorrection(quad.points, quad.points, s, ...
    soundSpeed, quad.weights, quad.outwardNormals());
K = K + Bw.' * (correction * quad.basis);
end


function rows = cqDoubleLayerPotential(surface, points, s, soundSpeed, quadratureOrder)
signs = surface.orientation.triangleOrientationSignsToOutward(:);
if any(signs == 0)
    error("volFemBemCoupledConvolutionQuadrature:orientation", ...
        "Surface orientation is unknown for %d triangle(s); cannot evaluate D(s).", ...
        sum(signs == 0));
end

tri = surface.tri;
vtx = surface.vtx;
rows = complex(zeros(size(points, 1), size(vtx, 1)));
for t = 1:size(tri, 1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t, :), :), points);
    rows(:, tri(t, :)) = rows(:, tri(t, :)) + signs(t) * J1 / (4 * pi);
end

quad = SurfaceQuadrature(surface, quadratureOrder);
correction = cqDoubleLayerCorrection(points, quad.points, s, ...
    soundSpeed, quad.weights, quad.outwardNormals());
rows = rows + correction * quad.basis;
end


function V = cqSingleLayerGalerkin(surface, s, soundSpeed, quadratureOrder)
quad = SurfaceQuadrature(surface, quadratureOrder);
nGauss = quad.nPoints();
nNodes = size(surface.vtx, 1);
tri = surface.tri;
vtx = surface.vtx;

P = zeros(nGauss, nNodes);
for t = 1:size(tri, 1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t, :), :), quad.points);
    P(:, tri(t, :)) = P(:, tri(t, :)) + I1;
end
Bw = quad.weightedBasis();
V = Bw.' * P / (4 * pi);

correction = cqSingleLayerCorrection(quad.points, quad.points, s, ...
    soundSpeed, quad.weights);
V = V + Bw.' * (correction * quad.basis);
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


function C = cqDoubleLayerCorrection(targetPoints, sourcePoints, s, soundSpeed, ...
        sourceWeights, sourceNormals)
nTarget = size(targetPoints, 1);
nSource = size(sourcePoints, 1);
C = complex(zeros(nTarget, nSource));
alpha = s / soundSpeed;
for i = 1:nTarget
    for j = 1:nSource
        delta = targetPoints(i, :) - sourcePoints(j, :);
        r = norm(delta);
        if r == 0
            value = 0.0;
        else
            normalDot = dot(delta, sourceNormals(j, :));
            base = normalDot / r^3;
            z = -alpha * r;
            value = base * stableExpTimesOneMinusZMinusOne(z);
        end
        C(i, j) = sourceWeights(j) * value / (4 * pi);
    end
end
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
            value = -alpha;
        else
            z = -alpha * r;
            value = stableExpm1OverR(z, r);
        end
        C(i, j) = sourceWeights(j) * value / (4 * pi);
    end
end
end


function value = stableExpTimesOneMinusZMinusOne(z)
%STABLEEXPTIMESONEMINUSZMINUSONE exp(z)*(1-z)-1, Taylor near z=0.

if abs(z) < 1e-5
    value = 0;
    for k = 2:10
        value = value + (1 - k) * z^k / factorial(k);
    end
else
    value = exp(z) * (1 - z) - 1;
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
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "four_tet_interior_node.vol"));
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
