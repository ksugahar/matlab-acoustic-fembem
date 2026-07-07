function tests = testAcousticRadiationForce
%TESTACOUSTICRADIATIONFORCE Radiation-force post-processor, series + BEM.
%
% The net acoustic radiation force on a rigid sphere in a +z plane wave,
% from the Brillouin radiation-stress control-surface integral. Gated the
% formula-free way (2026-07-04 measurements): the force is INDEPENDENT of
% the control radius to ~1e-10 (the divergence-free property - the primary
% correctness gate, no external formula needed), pushes downstream
% (F_z > 0), is axisymmetric (F_x = F_y = 0), and gives Y_p(kR=2) ~ 0.75.
% The SAME post-processor takes the analytic partial-wave series OR the BEM
% total field; BEM vs series agree to ~5% (faceting, the standard BEM
% band). This is the acoustic radiation-force / thrust physics that the
% adjoint AD will differentiate w.r.t. the phased-array phases.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testAnalyticSeriesForceIsControlRadiusIndependent(testCase)
k = 2.0;
pf = @(X) rigidSeriesField(k, X);
res = acousticRadiationForce(pf, k, "ControlRadius", 1.5);
verifyEqual(testCase, res.status, "ok");
verifyLessThan(testCase, res.controlRadiusResidual, 1e-6);   % measured 3e-10
verifyGreaterThan(testCase, res.force(3), 0);                % downstream
verifyLessThan(testCase, max(abs(res.force(1:2))), 1e-6 * res.force(3));
verifyEqual(testCase, res.forceFunction, 0.7515, "AbsTol", 0.02);
end


function testBemFieldForceMatchesSeries(testCase)
% the geometry-general payoff: the radiation force from the BEM-solved
% total field agrees with the analytic series to faceting (~5%), and is
% itself control-radius independent.
repoRoot = fileparts(fileparts(mfilename("fullpath")));
mesh = VolMesh(fullfile(repoRoot, "fixtures", "mesh_topology", ...
    "unit_sphere_fine.vol"));
surface = mesh.boundary();
k = 2.0;
sol = rigidScatteringSolve(surface, "Wavenumber", k, "QuadratureOrder", 7);
t = sol.trace;
pf = @(X) exp(1i*k*X(:,3)) + doubleLayerPotentialMatrix(surface, X, k, 7) * t;

res = acousticRadiationForce(pf, k, "ControlRadius", 1.5);
verifyLessThan(testCase, res.controlRadiusResidual, 1e-6);
verifyGreaterThan(testCase, res.force(3), 0);
verifyEqual(testCase, res.forceFunction, 0.7515, "AbsTol", 0.06);   % meas 0.715
end


function testPhysicalScaleAir40kHz(testCase)
% air @ 40 kHz, Rs = 2.73 mm sphere: the dimensionless Y_p converts to a
% physical force. Sanity: 160 dB SPL gives a sub-mN force (levitation /
% manipulation regime), 120 dB a sub-uN force.
k = 2.0;
res = acousticRadiationForce(@(X) rigidSeriesField(k, X), k, ...
    "ControlRadius", 1.5);
Y = res.forceFunction;

cAir = 343; rhoAir = 1.2; f = 40e3;
Rs = k / (2*pi*f/cAir);                 % 2.73 mm
verifyEqual(testCase, Rs, 2.73e-3, "AbsTol", 0.05e-3);
force = @(spl) Y * pi * Rs^2 * (20e-6 * 10^(spl/20))^2 / (rhoAir * cAir^2);
verifyLessThan(testCase, force(120), 1e-6);        % sub-uN
verifyGreaterThan(testCase, force(160), 1e-4);     % >0.1 mN
verifyLessThan(testCase, force(160), 1e-3);        % <1 mN
end


function p = rigidSeriesField(k, X)
% total field (incident + scattered) of a rigid unit sphere, +z plane wave
p = zeros(size(X, 1), 1);
for m = 1:size(X, 1)
    x = X(m, :);
    r = norm(x);
    ct = x(3) / r;
    L = ceil(k * r) + 15;
    Pp = 1; Pc = ct; sc = 0;
    for l = 0:L
        if l == 0
            Pl = Pp;
        elseif l == 1
            Pl = Pc;
        else
            Pl = ((2*l-1)*ct*Pc - (l-1)*Pp) / l; Pp = Pc; Pc = Pl;
        end
        dj = sphj(l-1, k) - (l+1)/k * sphj(l, k);
        dh = sphh(l-1, k) - (l+1)/k * sphh(l, k);
        sc = sc - (1i^l)*(2*l+1)*dj/dh * sphh(l, k*r) * Pl;
    end
    p(m) = exp(1i*k*x(3)) + sc;
end
end

function j = sphj(l, x), j = sqrt(pi/(2*x)) * besselj(l+0.5, x); end
function h = sphh(l, x), h = sqrt(pi/(2*x)) * (besselj(l+0.5, x) + 1i*bessely(l+0.5, x)); end
