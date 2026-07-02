function report = geometricIntegratorEnergyReport(stepSize, steps, omega)
%geometricIntegratorEnergyReport Harmonic-oscillator energy gate.
%
% This readable teaching helper compares three first lessons in time
% integration for q'' + omega^2 q = 0:
%
%   explicit Euler      : negative control, energy drifts upward
%   symplectic Euler    : geometric method, bounded energy error
%   implicit midpoint   : symplectic method, quadratic energy preserved here
%
% It returns only small structs and arrays so notebooks can jsonencode the
% result and the public radia gate can replay the same contract without MATLAB.

if nargin < 1 || isempty(stepSize)
    stepSize = 0.02;
end
if nargin < 2 || isempty(steps)
    steps = 1000;
end
if nargin < 3 || isempty(omega)
    omega = 1.0;
end

if ~(isscalar(stepSize) && isfinite(stepSize) && stepSize > 0)
    error("geometricIntegratorEnergyReport:stepSize", ...
        "stepSize must be a positive scalar.");
end
if ~(isscalar(steps) && isfinite(steps) && steps == round(steps) && steps > 0)
    error("geometricIntegratorEnergyReport:steps", ...
        "steps must be a positive integer.");
end
if ~(isscalar(omega) && isfinite(omega) && omega > 0)
    error("geometricIntegratorEnergyReport:omega", ...
        "omega must be a positive scalar.");
end

steps = double(steps);
q0 = 1.0;
p0 = 0.0;

explicit = integrateExplicitEuler(q0, p0, stepSize, steps, omega);
symplectic = integrateSymplecticEuler(q0, p0, stepSize, steps, omega);
midpoint = integrateImplicitMidpoint(q0, p0, stepSize, steps, omega);

rowExplicit = summarizeMethod("explicit_euler", explicit, stepSize, steps, omega);
rowSymplectic = summarizeMethod("symplectic_euler", symplectic, stepSize, steps, omega);
rowMidpoint = summarizeMethod("implicit_midpoint", midpoint, stepSize, steps, omega);
rows = [rowExplicit, rowSymplectic, rowMidpoint];

maxGeometricDrift = max([rowSymplectic.max_rel_energy_drift, rowMidpoint.max_rel_energy_drift]);
explicitToGeometricRatio = rowExplicit.max_rel_energy_drift / max(maxGeometricDrift, eps);

report = struct();
report.kind = "geometric_integrator_energy_report";
report.policy = "readable_geometric_time_integration_energy_gate";
report.problem = struct( ...
    "equation", "q'' + omega^2 q = 0", ...
    "hamiltonian", "0.5*(p^2 + omega^2*q^2)", ...
    "q0", q0, ...
    "p0", p0);
report.step_size_s = stepSize;
report.steps = steps;
report.omega_rad_per_s = omega;
report.method_rows = rows;
report.explicit_euler = rowExplicit;
report.symplectic_euler = rowSymplectic;
report.implicit_midpoint = rowMidpoint;
report.max_geometric_rel_energy_drift = maxGeometricDrift;
report.explicit_to_geometric_drift_ratio = explicitToGeometricRatio;
report.checks = struct( ...
    "symplecticEulerBounded", rowSymplectic.max_rel_energy_drift < 0.02, ...
    "implicitMidpointPreservesQuadraticEnergy", rowMidpoint.max_rel_energy_drift < 1e-10, ...
    "explicitEulerIsNegativeControl", explicitToGeometricRatio > 20);
report.pass = report.checks.symplecticEulerBounded ...
    && report.checks.implicitMidpointPreservesQuadraticEnergy ...
    && report.checks.explicitEulerIsNegativeControl;
end


function state = integrateExplicitEuler(q0, p0, h, steps, omega)
q = zeros(steps + 1, 1);
p = zeros(steps + 1, 1);
q(1) = q0;
p(1) = p0;
for k = 1:steps
    qOld = q(k);
    pOld = p(k);
    q(k + 1) = qOld + h * pOld;
    p(k + 1) = pOld - h * omega^2 * qOld;
end
state = struct("q", q, "p", p);
end


function state = integrateSymplecticEuler(q0, p0, h, steps, omega)
q = zeros(steps + 1, 1);
p = zeros(steps + 1, 1);
q(1) = q0;
p(1) = p0;
for k = 1:steps
    p(k + 1) = p(k) - h * omega^2 * q(k);
    q(k + 1) = q(k) + h * p(k + 1);
end
state = struct("q", q, "p", p);
end


function state = integrateImplicitMidpoint(q0, p0, h, steps, omega)
q = zeros(steps + 1, 1);
p = zeros(steps + 1, 1);
q(1) = q0;
p(1) = p0;
update = [1, -0.5 * h; 0.5 * h * omega^2, 1];
for k = 1:steps
    rhs = [q(k) + 0.5 * h * p(k); p(k) - 0.5 * h * omega^2 * q(k)];
    next = update \ rhs;
    q(k + 1) = next(1);
    p(k + 1) = next(2);
end
state = struct("q", q, "p", p);
end


function row = summarizeMethod(method, state, h, steps, omega)
energy = 0.5 * (state.p.^2 + (omega * state.q).^2);
energyInitial = energy(1);
energyFinal = energy(end);
relDrift = abs(energy - energyInitial) ./ energyInitial;

row = struct();
row.method = method;
row.steps = steps;
row.step_size_s = h;
row.omega_rad_per_s = omega;
row.energy_initial = energyInitial;
row.energy_final = energyFinal;
row.final_rel_energy_drift = abs(energyFinal - energyInitial) / energyInitial;
row.max_rel_energy_drift = max(relDrift);
end
