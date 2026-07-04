function res = elasticThrustAdjoint(model, sources, wavenumber, amplitudes, options)
%ELASTICTHRUSTADJOINT Radiation force on an ELASTIC bead and its sensitivity to
% phased-array amplitudes - wavefront-synthesised thrust through the FSI solve.
%
%   res = elasticThrustAdjoint(model, sources, k, amplitudes);
%   res.force            % 1x3 net radiation force on the elastic bead
%   res.objective        % the chosen thrust component (default F_z)
%   res.gradientReal     % dF/dRe(p_j), res.gradientImag = dF/dIm(p_j)
%   res.ascentDirection  % steepest-ascent complex direction = 2 Q p
%   res.forceForm        % {Q1,Q2,Q3}: reuse for a design loop, F_i = p^H Q_i p
%
% Closes the ultrasonic-thrust story: design the phased-array source amplitudes
% p (complex, one per source) to steer the acoustic RADIATION FORCE on a solid
% ELASTIC bead. Unlike a rigid/soft bead, the elastic interior has internal
% compressional + shear resonances (the FSI physics), so the force is computed
% through the full fluid-structure coupled solve (fsiCoupledSolve).
%
% Three ideas, each one readable step:
%
% 1. THE FIELD IS LINEAR IN THE AMPLITUDES. Solve the FSI once per source
%    (incident = a unit monopole at that source), keep its scattered surface
%    data. The total field at points X is then, for ALL sources at once,
%       fieldAt(X) = INC(X) + D_k(X) PS - S_k(X) QS        (nPts x nSrc)
%    where PS/QS are the per-source scattered pressure/flux on Gamma and D_k/S_k
%    are the double/single-layer potentials - two potential builds, SHARED over
%    sources (the incident does not change the surface operators). p_total = P p.
%
% 2. THE FORCE IS A QUADRATIC FORM. The radiation force is a SECOND-order
%    functional of the linear field, so each component is F_i(p) = p^H Q_i p
%    with Q_i Hermitian (F_i real). The field is SAMPLED once on the control
%    sphere; Q_i is assembled from the Brillouin radiation-stress flux.
%
% 3. THE GRADIENT IS ONE MATRIX-VECTOR PRODUCT. The Wirtinger steepest-ascent
%    direction of the real force is dF_i/dconj(p) = 2 Q_i p. So force AND its
%    full gradient are closed forms in p once Q_i is built - the N FSI solves
%    are the only real cost, and a design loop REUSES res.forceForm (no re-solve;
%    F_i(p) = real(p' * forceForm{i} * p), grad = 2 forceForm{i} p).
%
% Validated (Validate=true): the quadratic form reproduces a vectorised direct
% Brillouin integral (directForce) on the SAME field samples to machine
% precision at the design p AND a random p (independent assembly - the form is a
% Hermitian outer product, the direct integral forms the scalar field first);
% that direct integral matches the golden point-by-point acousticRadiationForce;
% the force is control-radius independent (div T = 0); and with GradientCheck the
% Wirtinger gradient matches central finite differences (exact for a quadratic).

arguments
    model (1,1) FemBemModel
    sources (:,3) double
    wavenumber (1,1) double {mustBePositive}
    amplitudes (:,1) double
    options.LongitudinalSpeed (1,1) double {mustBePositive} = 1.6
    options.ShearSpeed (1,1) double {mustBeNonnegative} = 0.9
    options.DensityRatio (1,1) double {mustBePositive} = 1.15
    options.FluidDensity (1,1) double {mustBePositive} = 1.0
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.ExteriorMethod (1,1) string ...
        {mustBeMember(options.ExteriorMethod, ["bem", "dtn"])} = "dtn"
    options.ForceComponent (1,1) double {mustBeMember(options.ForceComponent, [1 2 3])} = 3
    options.ControlRadius (1,1) double {mustBePositive} = 1.5
    options.NMu (1,1) double {mustBeInteger, mustBePositive} = 12
    options.NPhi (1,1) double {mustBeInteger, mustBePositive} = 24
    options.FieldQuadratureOrder (1,1) double {mustBeMember(options.FieldQuadratureOrder, [1 3 7])} = 3
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.AmplitudeStep (1,1) double {mustBePositive} = 1e-3
    options.Validate (1,1) logical = true
    options.GradientCheck (1,1) logical = false
end

k = wavenumber;
rho = options.FluidDensity;
c = options.SoundSpeed;
omega = k * c;
nSrc = size(sources, 1);
if numel(amplitudes) ~= nSrc
    error("elasticThrustAdjoint:amplitudes", ...
        "amplitudes must have one entry per source (%d).", nSrc);
end
p = amplitudes(:);
surface = model.surface;
Rc = options.ControlRadius;
h  = options.FiniteDifferenceStep;

% --- step 1: one FSI solve per source -> scattered surface data + incident ---
matArgs = {"Wavenumber", k, "LongitudinalSpeed", options.LongitudinalSpeed, ...
    "ShearSpeed", options.ShearSpeed, "DensityRatio", options.DensityRatio, ...
    "FluidDensity", rho, "ExteriorMethod", options.ExteriorMethod};
nB = size(surface.vtx, 1);
scatteredP = zeros(nB, nSrc);          % PS: per-source scattered pressure on Gamma
scatteredQ = zeros(nB, nSrc);          % QS: per-source scattered flux on Gamma
incidents = cell(1, nSrc);             % per-source incident-value handles
solStatus = strings(1, nSrc);
for j = 1:nSrc
    inc = monopoleIncident(sources(j, :), k);
    solj = fsiCoupledSolve(model, matArgs{:}, "Incident", inc);
    scatteredP(:, j) = solj.surfacePressure;
    scatteredQ(:, j) = solj.surfaceFlux;
    incidents{j} = inc.value;
    solStatus(j) = string(solj.status);
end

% the affine field map: p_total(X) = fieldAt(X) * p (two shared potential builds)
qOrder = options.FieldQuadratureOrder;
fieldAt = @(X) incidentValues(incidents, X) ...
    + doubleLayerPotentialMatrix(surface, X, k, qOrder) * scatteredP ...
    - singleLayerPotentialMatrix(surface, X, k, qOrder) * scatteredQ;

% --- step 2: sample the field on the control sphere, assemble the forms ---
sample = sampleFieldOnSphere(fieldAt, Rc, options.NMu, options.NPhi, h, omega, rho);
Qc = brillouinForms(sample, nSrc, rho, c);
force = [real(p' * Qc{1} * p), real(p' * Qc{2} * p), real(p' * Qc{3} * p)];

% --- step 3: the chosen thrust component and its Wirtinger gradient ---
ci = options.ForceComponent;
Q = Qc{ci};
objective = real(p' * Q * p);
gradientReal = 2 * real(Q * p);
gradientImag = 2 * imag(Q * p);
ascentDirection = 2 * (Q * p);                 % dF/dconj(p), steepest ascent

res = struct();
res.kind = "elastic_bead_radiation_force_phased_array_adjoint";
res.policy = "quadratic_form_of_the_fsi_linear_field_wirtinger_gradient";
res.wavenumber = k;
res.numSources = nSrc;
res.fsiSolves = nSrc;
res.exteriorMethod = options.ExteriorMethod;
res.forceComponent = ci;
res.force = force;
res.objective = objective;
res.gradientReal = gradientReal;
res.gradientImag = gradientImag;
res.ascentDirection = ascentDirection;
res.forceForm = Qc;                            % reuse in a design loop (no re-solve)
res.solveStatus = solStatus;

if ~options.Validate
    res.status = ternary(all(solStatus == "ok"), "ok", "needs_attention");
    return
end

% --- validation ---
% (a) the quadratic form vs a vectorised direct integral on the SAME samples
%     (independent assembly: outer-product form vs scalar-field integral)
consErr = 0;
for a = {p, probeAmplitudes(nSrc, 1)}
    fQ = [real(a{1}' * Qc{1} * a{1}), real(a{1}' * Qc{2} * a{1}), real(a{1}' * Qc{3} * a{1})];
    fD = directForce(sample, a{1}, rho, c);
    consErr = max(consErr, norm(fQ - fD) / max(norm(fD), eps));
end
res.consistencyError = consErr;

% (b) the direct integral vs the golden point-by-point acousticRadiationForce
%     (a tiny shared quadrature - only checks the two integrators agree)
tiny = sampleFieldOnSphere(fieldAt, Rc, 4, 6, h, omega, rho);
fRef = acousticRadiationForce(@(X) fieldAt(X) * p, k, "ControlRadius", Rc, ...
    "NMu", 4, "NPhi", 6, "FiniteDifferenceStep", h, "Rho", rho, "SoundSpeed", c).force;
res.independentForceError = norm(directForce(tiny, p, rho, c) - fRef) / max(norm(fRef), eps);

% (c) control-radius independence (div T = 0): force at Rc vs 1.4 Rc
outer = sampleFieldOnSphere(fieldAt, 1.4 * Rc, options.NMu, options.NPhi, h, omega, rho);
res.controlRadiusResidual = norm(force - directForce(outer, p, rho, c)) / max(norm(force), eps);

checks = struct( ...
    "fsiSolvesOk", all(solStatus == "ok"), ...
    "formMatchesDirect", consErr < 1e-8, ...
    "directMatchesGolden", res.independentForceError < 1e-8, ...
    "controlRadiusIndependent", res.controlRadiusResidual < 5e-3);   % quadrature-limited (3e-6 at 12/24)

% (d) Wirtinger gradient vs central FD of the (validated) force form
if options.GradientCheck
    hg = options.AmplitudeStep;
    forceComp = @(a) real(a' * Q * a);
    fdRe = zeros(nSrc, 1);
    fdIm = zeros(nSrc, 1);
    for j = 1:nSrc
        e = zeros(nSrc, 1); e(j) = 1;
        fdRe(j) = (forceComp(p + hg * e) - forceComp(p - hg * e)) / (2 * hg);
        fdIm(j) = (forceComp(p + 1i * hg * e) - forceComp(p - 1i * hg * e)) / (2 * hg);
    end
    res.gradientCheckRelError = max( ...
        norm(gradientReal - fdRe) / max(norm(fdRe), eps), ...
        norm(gradientImag - fdIm) / max(norm(fdIm), eps));
    checks.gradientMatchesFiniteDifference = res.gradientCheckRelError < 1e-6;
end

res.checks = checks;
res.status = ternary(all(structfun(@(v) logical(v), checks)), "ok", "needs_attention");
end


% =====================================================================
function s = sampleFieldOnSphere(fieldAt, Rc, nMu, nPhi, h, omega, rho)
%SAMPLEFIELDONSPHERE Amplitude-sensitivities of the pressure and velocity at the
% control-sphere quadrature points (the only place fieldAt is evaluated in bulk:
% the point set and its six one-sided FD neighbours). Returns a struct with
%   pSens : nQ x nSrc   d p_total / d p
%   vSens : 1x3 cell    d v_i / d p  (velocity = grad p / (i omega rho), Euler)
%   normal, weight, point : the Gauss-Legendre x phi control-sphere quadrature.
[cosTheta, wTheta] = gaussLegendreNodes(nMu);
phi  = (0:nPhi - 1).' / nPhi * 2 * pi;
dPhi = 2 * pi / nPhi;

[iTheta, iPhi] = ndgrid(1:nMu, 1:nPhi);
iTheta = iTheta(:); iPhi = iPhi(:);
sinTheta = sqrt(1 - cosTheta(iTheta).^2);
s.normal = [sinTheta .* cos(phi(iPhi)), sinTheta .* sin(phi(iPhi)), cosTheta(iTheta)];
s.point  = Rc * s.normal;
s.weight = wTheta(iTheta) * dPhi * Rc^2;

s.pSens = fieldAt(s.point);
s.vSens = cell(1, 3);
for i = 1:3
    step = zeros(1, 3); step(i) = h;
    s.vSens{i} = (fieldAt(s.point + step) - fieldAt(s.point - step)) / (2 * h) / (1i * omega * rho);
end
end


function Qc = brillouinForms(s, nSrc, rho, c)
%BRILLOUINFORMS Radiation-force components as Hermitian quadratic forms in the
% amplitudes, F_i(p) = p^H Q_i p, from the field samples. The pressure/velocity
% are linear in p, so the Brillouin stress flux is quadratic in p.
%
% Sensitivity convention: p_total = dp.' * p, so |p_total|^2 = p^H (conj(dp) dp.')
% p - the Hermitian outer product is conj(z) z.', NOT z z' (they differ for
% COMPLEX amplitudes; z z' would silently give |dp' p|^2).
outerForm = @(z) conj(z) * z.';                     % |z.' p|^2 = p^H outerForm(z) p

Q = {zeros(nSrc), zeros(nSrc), zeros(nSrc)};
for q = 1:size(s.point, 1)
    dp = s.pSens(q, :).';
    dv = {s.vSens{1}(q, :).', s.vSens{2}(q, :).', s.vSens{3}(q, :).'};
    n  = s.normal(q, :);

    % Lagrangian density L = |p|^2/(4 rho c^2) - rho |v|^2 / 4 (a Hermitian form)
    lagrangian = outerForm(dp) / (4 * rho * c^2) ...
        - (rho / 4) * (outerForm(dv{1}) + outerForm(dv{2}) + outerForm(dv{3}));
    normalVel = dv{1} * n(1) + dv{2} * n(2) + dv{3} * n(3);   % (v.n) sensitivity

    for i = 1:3
        % (T n)_i = L n_i + (rho/2) Re[ v_i conj(v.n) ], each a Hermitian form
        reynolds  = conj(normalVel) * dv{i}.';
        integrand = n(i) * lagrangian + (rho / 4) * (reynolds + reynolds');
        Q{i} = Q{i} - s.weight(q) * integrand;
    end
end

Qc = {(Q{1} + Q{1}') / 2, (Q{2} + Q{2}') / 2, (Q{3} + Q{3}') / 2};   % Hermitian
end


function F = directForce(s, p, rho, c)
%DIRECTFORCE Brillouin radiation force for amplitudes p from the field samples -
% the SCALAR-field integral (form the field first, then the stress), an
% assembly independent of brillouinForms' outer products (so it catches a form
% conjugation error). Vectorised twin of acousticRadiationForce.
pressure = s.pSens * p;                             % nQ x 1 total pressure
vel = [s.vSens{1} * p, s.vSens{2} * p, s.vSens{3} * p];   % nQ x 3 velocity
F = [0 0 0];
for q = 1:size(s.point, 1)
    v = vel(q, :);
    lagrangian = abs(pressure(q))^2 / (4 * rho * c^2) - rho * (v * v') / 4;
    T = real(lagrangian) * eye(3) + (rho / 2) * real(v.' * conj(v));
    F = F - (T * s.normal(q, :).').' * s.weight(q);
end
end


% =====================================================================
function V = incidentValues(incidents, X)
%INCIDENTVALUES Per-source incident pressures as columns [f_1(X) ... f_n(X)].
V = cell2mat(cellfun(@(f) f(X), incidents, "UniformOutput", false));
end


function inc = monopoleIncident(src, k)
%MONOPOLEINCIDENT Incident struct (value + gradient) of a point source at src.
inc = struct("value", @(X) monopoleValue(X, src, k), "grad", @(X) monopoleGrad(X, src, k));
end

function p = monopoleValue(X, src, k)
r = sqrt(sum((X - src).^2, 2));
p = exp(1i * k * r) ./ (4 * pi * r);
end

function g = monopoleGrad(X, src, k)
d = X - src;
r = sqrt(sum(d.^2, 2));
pp = exp(1i * k * r) ./ (4 * pi * r);
g = (pp .* (1i * k - 1 ./ r) ./ r) .* d;         % grad = p (ik - 1/r) d/r
end


function v = probeAmplitudes(n, seed)
%PROBEAMPLITUDES Deterministic complex amplitudes to pin the force form.
rng(seed);
v = randn(n, 1) + 1i * randn(n, 1);
end


function [x, w] = gaussLegendreNodes(n)
%GAUSSLEGENDRENODES Golub-Welsch nodes/weights on [-1, 1] (matches
% acousticRadiationForce's control-sphere quadrature).
beta = (1:n-1) ./ sqrt(4 * (1:n-1).^2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[x, ix] = sort(diag(D));
w = 2 * (V(1, ix).^2).';
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
