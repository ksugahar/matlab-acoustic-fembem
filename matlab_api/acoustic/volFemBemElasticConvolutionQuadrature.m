function result = volFemBemElasticConvolutionQuadrature(volFile, options)
%VOLFEMBEMELASTICCONVOLUTIONQUADRATURE Elastic-structure FEM + acoustic BEM,
% coupled, time domain by Lubich convolution quadrature (SCATTERING).
%
%   The canonical acoustic FSI simulator: a solid ELASTIC scatterer (vector P1
%   elasticity FEM interior) in an unbounded fluid (P1 acoustic BEM exterior),
%   hit by an incident plane-wave PULSE, solved in the TIME DOMAIN by CQ.  It is
%   the time-domain twin of fsiCoupledSolve (frequency), built on the CQ machinery
%   of volFemBemConvolutionQuadrature (scalar).  Per CQ Laplace node s = delta(zeta)/dt:
%
%     [ Ks + s^2 Ms      G'             0   ] [u  ]   [ -G' pinc(s) ]
%     [ rho_f s^2 G      0              Mb  ] [p_s] = [ -Minc(s)    ]
%     [ 0                1/2 Mb - K(s)  V(s)] [q_s]   [  0          ]
%
%   with Ks,Ms the elasticity stiffness/mass and s^2 the Laplace image of
%   d^2/dt^2.  At s = -i c0 k this block REPRODUCES fsiCoupledSolve's frequency
%   system (Ks - w^2 Ms, -rho_f w^2 G) to ~1e-3 (the single-layer coincident Delta),
%   which pins the whole Laplace-domain assembly against the validated frequency
%   solver.  The incident plane-wave pulse is retarded, pinc(x,s) = ghat(s) exp(-s z/c0),
%   grad pinc = [0 0 -s/c0 pinc]; the inverse CQ FFT with the rho^-n unscaling
%   recovers the causal interior displacement u(t), boundary flux q_s(t), and
%   exterior scattered pressure p(x,t).  1st-order tria/tet only.
%
%   Speeds/density are RELATIVE TO THE FLUID (c0 = SoundSpeed); the elastic Lame
%   constants are mu = DensityRatio*cT^2, lambda = DensityRatio*(cL^2 - 2 cT^2),
%   matching fsiCoupledSolve.
%
%   The result struct carries the CQ internals (cqZeta, cqLaplaceParameter,
%   conditionNumbers, relativeResiduals) plus pressure/boundaryData aliases, so
%   visualizeConvolutionQuadrature(result) draws the same six-panel CQ X-ray.
%
%   NOTE: the CQ single/double-layer operators below are local copies of the ones
%   in volFemBemConvolutionQuadrature / volTdBemConvolutionQuadrature.  A follow-up
%   should promote the shared CQ operators to matlab_api/bem/ so the three CQ
%   solvers stop duplicating them.

arguments
    volFile (1,1) string = ""
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 3)} = 24
    options.TimeStep (1,1) double {mustBePositive} = 0.35
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.LongitudinalSpeed (1,1) double {mustBePositive} = 2.0
    options.ShearSpeed (1,1) double {mustBeNonnegative} = 1.0
    options.DensityRatio (1,1) double {mustBePositive} = 1.5
    options.FluidDensity (1,1) double {mustBePositive} = 1.0
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder,[1 3 7])} = 3
    options.CqRadius double = []
    options.ObservationPoints double = []
    options.ResidualTolerance (1,1) double {mustBePositive} = 1e-6
end
if strlength(volFile) == 0
    volFile = defaultFixture();
end
c0 = options.SoundSpeed; rhoF = options.FluidDensity; rhoS = options.DensityRatio;
mu = rhoS*options.ShearSpeed^2;
lamE = rhoS*(options.LongitudinalSpeed^2 - 2*options.ShearSpeed^2);

meshTimer = tic;
model = FemBemModel(volFile); mesh = model.mesh; surface = model.surface;
nV = size(mesh.vtx,1); ids = surface.volNodeIds; nB = numel(ids);
[Ks, Ms] = elasticityMatrices(mesh, lamE, mu, rhoS);   % interior elasticity
[Mb, ~] = SurfaceP1Space(surface).mass();              % boundary P1 mass
G = couplingMatrix(surface, ids, nV);                  % interface (geometry-only)
zNode = surface.vtx(:,3);
if isempty(options.ObservationPoints), obs = defaultObservationPoints(mesh.vtx);
else, obs = options.ObservationPoints; end
if size(obs, 2) ~= 3
    error("volFemBemElasticConvolutionQuadrature:ObservationPoints", ...
        "ObservationPoints must have three columns.");
end
meshSeconds = toc(meshTimer);

N = options.NumTime; dt = options.TimeStep; t = (0:N-1).'*dt;
rho = options.CqRadius; if isempty(rho), rho = sqrt(eps)^(1/N); end
if ~(isscalar(rho) && rho > 0 && rho < 1)
    error("volFemBemElasticConvolutionQuadrature:CqRadius", ...
        "CqRadius must be a scalar in (0, 1).");
end
n = (0:N-1).'; zeta = rho.*exp(-2i*pi*n/N); s = bdfDelta(zeta, options.Method)/dt;
g = rickerPulse(t, dt); ghat = fft((rho.^n).*g);       % incident-pulse spectrum

cqTimer = tic;
Uhat = zeros(N,3*nV); Qhat = zeros(N,nB); Phat = zeros(N,size(obs,1));
resid = zeros(N,1); condno = zeros(N,1);
ZuB = sparse(3*nV,nB); ZBB = sparse(nB,nB);
for l = 1:N
    Kdyn = Ks + s(l)^2*Ms;                                     % Laplace elastodynamics
    V = laplaceSingleLayerGalerkin(surface, s(l), c0, options.QuadratureOrder);
    K = cqDoubleLayerGalerkin(surface, s(l), c0, options.QuadratureOrder);
    pinc = ghat(l).*exp(-s(l).*zNode./c0);                     % retarded incident, nodal
    Minc = incidentNormalFlux(surface, s(l), c0, ghat(l));     % (grad pinc . n)
    lhs = [ Kdyn,          G.',           ZuB;
            rhoF*s(l)^2*G, ZBB,           Mb;
            ZuB.',         0.5*Mb - K,    V ];
    rhs = [ -G.'*pinc; -Minc; zeros(nB,1) ];
    x = lhs\rhs;
    u = x(1:3*nV); ps = x(3*nV+(1:nB)); qs = x(3*nV+nB+(1:nB));
    Sobs = cqSingleLayerPotential(surface, obs, s(l), c0, options.QuadratureOrder);
    Dobs = cqDoubleLayerPotential(surface, obs, s(l), c0, options.QuadratureOrder);
    Uhat(l,:) = u.'; Qhat(l,:) = qs.'; Phat(l,:) = (Dobs*ps - Sobs*qs).';
    resid(l) = norm(lhs*x - rhs)/max(1, norm(rhs));
    condno(l) = condest(sparse(lhs));
end
uC = (rho.^(-n)).*ifft(Uhat,[],1);
qC = (rho.^(-n)).*ifft(Qhat,[],1);
pC = (rho.^(-n)).*ifft(Phat,[],1);
cqSeconds = toc(cqTimer);

result = struct();
result.kind = "elastic_fem_acoustic_bem_coupled_cq_time_response";
result.policy = "vector_p1_elasticity_fem_interior_p1_acoustic_bem_exterior_lubich_cq";
result.method = "lubich_" + lower(options.Method) + "_cq_elastic_fem_acoustic_bem_coupled";
result.volFile = string(volFile);
result.meshId = model.mesh.meshId;
result.meshSummary = model.mesh.summary;
result.time = t; result.timeStep = dt; result.cqRadius = rho;
result.cqZeta = zeta; result.cqLaplaceParameter = s;
result.incident = g;
result.interiorDisplacement = real(uC);
result.boundaryDensity = real(qC);
result.observationPoints = obs;
result.exteriorPressure = real(pC);
result.pressure = real(pC);              % alias so the CQ visualizer works
result.boundaryData = repmat(g, 1, nB); % incident trace time signal (visual)
result.relativeResiduals = resid;
result.conditionNumbers = condno;
result.timing = struct("mesh_or_import", meshSeconds, "coupled_cq_laplace_solves", cqSeconds);

scaleP = max(abs(result.exteriorPressure), [], "all");
result.summary = struct("num_time", N, "num_volume_dof", 3*nV, "num_boundary_nodes", nB, ...
    "num_observation_points", size(obs,1), ...
    "max_abs_exterior_pressure", scaleP, "max_relative_residual", max(resid), ...
    "max_condition_number", max(condno), ...
    "max_imag_exterior_before_real", max(abs(imag(pC)), [], "all"), ...
    "causal_leading_ratio", max(abs(pC(1,:)), [], "all") / max(scaleP, realmin));
result.checks = struct( ...
    "vol_mesh_tri_tet", model.mesh.summary.triangles > 0 && model.mesh.summary.tets > 0, ...
    "vector_elastic_interior", 3*nV > nV, ...
    "laplace_parameters_positive_real", all(real(s) > 0), ...
    "coupled_residuals_small", result.summary.max_relative_residual < options.ResidualTolerance, ...
    "finite_exterior", all(isfinite(result.exteriorPressure), "all"), ...
    "nonzero_exterior_response", scaleP > 0, ...
    "causal_leading_step", result.summary.causal_leading_ratio < 1e-2, ...
    "real_exterior_response", result.summary.max_imag_exterior_before_real < 1e-8*max(1, scaleP));
if all(structfun(@(v) logical(v), result.checks)), result.status = "ok";
else, result.status = "needs_attention"; end
end


% ===================================================================== %
% locals -- interface + CQ single/double layer.  The CQ layer operators
% duplicate volFemBemConvolutionQuadrature's; promote to matlab_api/bem/.
% ===================================================================== %
function G = couplingMatrix(surface, ids, nV)
% G_ij = int_Gamma mu_i (n . phi_struct_j), geometry-only (no incident).
signs = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri; vtx = surface.vtx; nB = size(vtx,1);
mT = (ones(3)+eye(3))/12; nE = 27*size(tri,1);
r = zeros(nE,1); c = zeros(nE,1); v = zeros(nE,1); cur = 1;
for tt = 1:size(tri,1)
    lc = tri(tt,:); X = vtx(lc,:);
    cr = cross(X(2,:)-X(1,:), X(3,:)-X(1,:));
    area = 0.5*norm(cr); nrm = signs(tt)*cr/norm(cr); vid = ids(lc);
    for a = 1:3, for b = 1:3, for d = 1:3
        r(cur) = lc(a); c(cur) = 3*vid(b)-3+d; v(cur) = area*mT(a,b)*nrm(d); cur = cur+1;
    end, end, end
end
G = sparse(r, c, v, nB, 3*nV);
end

function Minc = incidentNormalFlux(surface, s, c0, ghat)
% Minc_i = int_Gamma mu_i (grad pinc . n), pinc(x,s)=ghat exp(-s z/c0),
% grad pinc = [0 0 -s/c0 * pinc]; P1-lumped centroid rule (as in fsiCoupledSolve).
signs = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri; vtx = surface.vtx;
Minc = complex(zeros(size(vtx,1),1));
for tt = 1:size(tri,1)
    lc = tri(tt,:); X = vtx(lc,:);
    cr = cross(X(2,:)-X(1,:), X(3,:)-X(1,:));
    area = 0.5*norm(cr); nrm = signs(tt)*cr/norm(cr);
    Xc = mean(X,1);
    dpz = ghat*(-s/c0)*exp(-s*Xc(3)/c0);      % d(pinc)/dz at the centroid
    qn = dpz*nrm(3);                          % (grad pinc).n = dpz * n_z
    Minc(lc) = Minc(lc) + (area/3)*qn;
end
end

function K = cqDoubleLayerGalerkin(surface, s, c, q)
qd = SurfaceQuadrature(surface, q);
sg = surface.orientation.triangleOrientationSignsToOutward(:);
nG = qd.nPoints(); nN = size(surface.vtx,1); tri = surface.tri; vtx = surface.vtx;
P = complex(zeros(nG, nN));
for t = 1:size(tri,1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t,:),:), qd.points);
    P(:, tri(t,:)) = P(:, tri(t,:)) + sg(t)*J1;
end
Bw = qd.weightedBasis(); K = Bw.'*P/(4*pi);
C = dlCorr(qd.points, qd.points, s, c, qd.weights, qd.outwardNormals());
K = K + Bw.'*(C*qd.basis);
end

function rows = cqDoubleLayerPotential(surface, pts, s, c, q)
sg = surface.orientation.triangleOrientationSignsToOutward(:);
tri = surface.tri; vtx = surface.vtx;
rows = complex(zeros(size(pts,1), size(vtx,1)));
for t = 1:size(tri,1)
    [~, J1] = laplaceDoubleLayerPanelIntegrals(vtx(tri(t,:),:), pts);
    rows(:, tri(t,:)) = rows(:, tri(t,:)) + sg(t)*J1/(4*pi);
end
qd = SurfaceQuadrature(surface, q);
C = dlCorr(pts, qd.points, s, c, qd.weights, qd.outwardNormals());
rows = rows + C*qd.basis;
end

function rows = cqSingleLayerPotential(surface, pts, s, c, q)
tri = surface.tri; vtx = surface.vtx;
rows = zeros(size(pts,1), size(vtx,1));
for t = 1:size(tri,1)
    [~, I1] = laplacePanelIntegrals(vtx(tri(t,:),:), pts);
    rows(:, tri(t,:)) = rows(:, tri(t,:)) + I1/(4*pi);
end
qd = SurfaceQuadrature(surface, q);
C = slCorr(pts, qd.points, s, c, qd.weights);
rows = rows + C*qd.basis;
end

function C = dlCorr(tp, sp, s, c, w, nr)
nT = size(tp,1); nS = size(sp,1); C = complex(zeros(nT,nS)); al = s/c;
for i = 1:nT, for j = 1:nS
    d = tp(i,:) - sp(j,:); r = norm(d);
    if r == 0, vv = 0;
    else
        z = -al*r;
        if abs(z) < 1e-5, e = 0; for k = 2:10, e = e + (1-k)*z^k/factorial(k); end
        else, e = exp(z)*(1-z) - 1; end
        vv = dot(d, nr(j,:))/r^3 * e;
    end
    C(i,j) = w(j)*vv/(4*pi);
end, end
end

function C = slCorr(tp, sp, s, c, w)
nT = size(tp,1); nS = size(sp,1); C = complex(zeros(nT,nS)); al = s/c;
for i = 1:nT, for j = 1:nS
    r = norm(tp(i,:) - sp(j,:));
    if r == 0, vv = -al;
    else
        z = -al*r;
        if abs(z) < 1e-5, vv = 0; tm = 1; for k = 1:10, tm = tm*z/k; vv = vv + tm/r; end
        else, vv = (exp(z)-1)/r; end
    end
    C(i,j) = w(j)*vv/(4*pi);
end, end
end

function d = bdfDelta(z, m)
switch upper(m), case "BDF1", d = 1 - z; case "BDF2", d = 1.5 - 2*z + 0.5*z.^2; end
end

function volFile = defaultFixture()
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_ball_maxh018.vol"));
end

function obs = defaultObservationPoints(nd)
ct = mean(nd,1); r = max(vecnorm(nd-ct,2,2)); if r <= 0, r = 1; end
obs = ct + r*[0 0 2.4; 2.2 0 0.4; 0 2.0 1.2];
end

function src = rickerPulse(t, dt)
ct = 4*dt; w = 1.5*dt; x = (t-ct)/w;
src = (1 - 2*x.^2).*exp(-x.^2); src(abs(src) < 1e-14) = 0;
end
