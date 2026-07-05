function res = acousticRadiationForce(pressureField, wavenumber, options)
%ACOUSTICRADIATIONFORCE Net acoustic radiation force by the stress integral.
%
%   res = acousticRadiationForce(pressureField, k);
%   res.force              % 1x3 net radiation force on the scatterer
%   res.forceFunction      % Y_p = F_z / (pi Rs^2 <E>), the dimensionless
%                          %       radiation force function
%   res.controlRadiusResidual  % self-consistency: |F| must not depend on
%                              % where the control surface is placed
%
% The time-averaged acoustic radiation force is a SECOND-order functional of
% the LINEAR field: the closed-surface integral of the Brillouin
% radiation-stress tensor over any control sphere enclosing the scatterer,
%
%   T_ij = [ |p|^2/(4 rho c^2) - rho |v|^2 / 4 ] delta_ij
%          + (rho/2) Re[v_i conj(v_j)],       v = grad p / (i omega rho),
%   F_i  = -oint_S T_ij n_j dS   (n outward).
%
% T is divergence-free in the source-free fluid, so F is INDEPENDENT of the
% control radius - the primary, formula-free correctness gate (this
% function evaluates at ControlRadius and 1.4*ControlRadius and reports the
% residual). pressureField is any callable X (Npts x 3) -> complex pressure
% column, so the SAME post-processor takes the analytic partial-wave series
% OR a BEM total field (incident + D_k[trace]) - geometry-general.
%
% Dimensionless by default (rho = c = 1, omega = k, unit scatterer and
% incident amplitude); the physical force in a medium is
%   F_z = Y_p * pi * Rs^2 * <E>,   <E> = p_rms^2 / (rho c^2) = |A|^2 / (2 rho c^2)
%                                        for a progressive wave (p_rms = |A|/sqrt2).
% Air @ 40 kHz (c = 343, rho = 1.2, lambda = 8.58 mm): kR = 2 is a
% Rs = 2.73 mm scatterer; at 140 dB SPL (<E> = 0.28 J/m^3) a Y_p ~ 0.75
% sphere feels ~5 uN, at 160 dB ~0.5 mN - the acoustic-manipulation regime.
%
% Reference gate: for a rigid sphere in a plane progressive wave the force
% pushes downstream (F_z > 0) and Y_p(kR=2) ~ 0.75 (King 1934).

arguments
    pressureField (1,1) function_handle
    wavenumber (1,1) double {mustBePositive}
    options.ControlRadius (1,1) double {mustBePositive} = 1.5
    options.NMu (1,1) double {mustBeInteger, mustBePositive} = 20
    options.NPhi (1,1) double {mustBeInteger, mustBePositive} = 40
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.Rho (1,1) double {mustBePositive} = 1.0
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.ScattererRadius (1,1) double {mustBePositive} = 1.0
    options.IncidentAmplitude (1,1) double = 1.0
end

k = wavenumber;
rho = options.Rho;
c = options.SoundSpeed;
omega = k * c;

F1 = stressIntegral(pressureField, options.ControlRadius, k, rho, c, omega, options);
F2 = stressIntegral(pressureField, 1.4 * options.ControlRadius, k, rho, c, omega, options);

energyDensity = abs(options.IncidentAmplitude)^2 / (2 * rho * c^2);
Y = F1(3) / (pi * options.ScattererRadius^2 * energyDensity);

res = struct();
res.kind = "acoustic_radiation_force_brillouin_stress_integral";
res.policy = "second_order_functional_of_the_linear_field_control_surface";
res.wavenumber = k;
res.force = F1;
res.forceSecondRadius = F2;
res.controlRadiusResidual = norm(F1 - F2) / max(norm(F1), eps);
res.forceFunction = Y;
res.energyDensity = energyDensity;
res.controlRadius = options.ControlRadius;
res.checks = struct( ...
    "controlRadiusIndependent", res.controlRadiusResidual < 1e-4, ...
    "axisymmetric", max(abs(F1(1:2))) < 1e-2 * abs(F1(3)), ...
    "pushesDownstream", F1(3) > 0);
if all(structfun(@(v) logical(v), res.checks))
    res.status = "ok";
else
    res.status = "needs_attention";
end
end


function F = stressIntegral(pressureField, Rc, k, rho, c, omega, options)
%STRESSINTEGRAL Brillouin radiation-stress flux over a control sphere.
[mu, wmu] = gaussLegendreNodes(options.NMu);
phi = (0:options.NPhi - 1).' / options.NPhi * 2 * pi;
dphi = 2 * pi / options.NPhi;
h = options.FiniteDifferenceStep;

F = [0 0 0];
for a = 1:options.NMu
    st = sqrt(1 - mu(a)^2);
    for b = 1:options.NPhi
        n = [st * cos(phi(b)), st * sin(phi(b)), mu(a)];   % outward radial
        x = Rc * n;
        p = pressureField(x);
        gp = zeros(1, 3);
        for i = 1:3
            e = zeros(1, 3);
            e(i) = h;
            gp(i) = (pressureField(x + e) - pressureField(x - e)) / (2 * h);
        end
        v = gp / (1i * omega * rho);                       % Euler relation
        lagrangian = abs(p)^2 / (4 * rho * c^2) - rho * (v * v') / 4;
        T = real(lagrangian) * eye(3) + (rho / 2) * real(v.' * conj(v));
        F = F - (T * n.').' * (wmu(a) * dphi * Rc^2);
    end
end
end


function [x, w] = gaussLegendreNodes(n)
%GAUSSLEGENDRENODES Golub-Welsch nodes/weights on [-1, 1].
beta = (1:n-1) ./ sqrt(4 * (1:n-1).^2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[x, ix] = sort(diag(D));
w = 2 * (V(1, ix).^2).';
end
