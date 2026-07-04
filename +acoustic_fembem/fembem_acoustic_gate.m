function report = fembem_acoustic_gate(kind, options)
%FEMBEM_ACOUSTIC_GATE Run an integrated Gypsilab acoustic solve vs analytic.
%
%   report = acoustic_fembem.fembem_acoustic_gate("soft");
%   report = acoustic_fembem.fembem_acoustic_gate("rigid", Wavenumber=2.0);
%   report = acoustic_fembem.fembem_acoustic_gate("coupled_fluid");
%
% The FEM/BEM analogue of acoustic_fembem.result_manifest_gate: a runnable
% correctness gate that drives the integrated Gypsilab solver
% (acoustic_fembem.repository_root) on a committed sphere/ball fixture and compares the
% result to the analytic partial-wave series, returning a struct verdict.
% Pure MATLAB - no NGSolve, no committed .mat needed (those back the
% cross-code legs; this gate is the analytic leg).
%
% kind:
%   "soft"              exterior sound-soft scattering vs the soft-sphere
%                       series (singleLayerDirichletSolve)
%   "rigid"             sound-hard total-field solve vs the rigid-sphere
%                       series (rigidScatteringSolve)
%   "coupled_invisible" coupled transmission with k1=k0, rho1=rho0 - the
%                       exact invisibility null gate (interior == incident)
%   "coupled_fluid"     coupled transmission vs the Anderson fluid-sphere
%                       series (c1/c0 = 0.7, rho1/rho0 = 1.2)
%   "focus_adjoint"     reverse-mode adjoint AD: the phased-array wavefront-
%                       synthesis gradient vs central finite differences
%                       (relative_error = the adjoint-vs-FD gradient check)
%   "radiation_force"   acoustic radiation force on a rigid sphere from the
%                       BEM field vs the analytic series (ultrasonic thrust;
%                       relative_error = |Y_bem - Y_series| / Y_series)
%   "fsi"               fluid-structure interaction: the elastic-sphere FSI
%                       coupled solve (vector elasticity FEM + acoustic BEM)
%                       stiff limit vs the rigid sphere (formulation gate)

arguments
    kind (1,1) string {mustBeMember(kind, ...
        ["soft", "rigid", "coupled_invisible", "coupled_fluid", ...
         "focus_adjoint", "radiation_force", "fsi", "fsi_dtn", "thrust_adjoint"])}
    options.Fixture (1,1) string = ""
    options.Wavenumber (1,1) double {mustBePositive} = 2.0
    options.QuadratureOrder (1,1) double {mustBeMember( ...
        options.QuadratureOrder, [1 3 7])} = 7
    options.Tolerance (1,1) double = NaN
    options.Probes (:,3) double = [2 0 0; 0 0 3; -1.2 1.6 0]
end

root = acoustic_fembem.repository_root();
addpath(genpath(fullfile(root, "matlab_api")));

isCoupled = startsWith(kind, "coupled") || startsWith(kind, "fsi") || kind == "thrust_adjoint";
fixture = options.Fixture;
if fixture == ""
    if isCoupled
        fixture = "unit_ball_maxh018.vol";
    else
        fixture = "unit_sphere_fine.vol";
    end
end
volPath = fullfile(root, "fixtures", "mesh_topology", fixture);
if ~isfile(volPath)
    error("acoustic_fembem:fembemFixtureMissing", ...
        "Fixture not found: %s", volPath);
end

k = options.Wavenumber;
probes = options.Probes;
tol = options.Tolerance;
defaultTol = struct("soft", 0.06, "rigid", 0.08, ...
    "coupled_invisible", 0.07, "coupled_fluid", 0.10, ...
    "focus_adjoint", 1e-6, "radiation_force", 0.06, "fsi", 1e-2, ...
    "fsi_dtn", 1e-2, "thrust_adjoint", 1e-6);
if isnan(tol)
    tol = defaultTol.(kind);
end

report = struct();
report.tool = "acoustic_fembem_acoustic_gate";
report.kind = kind;
report.fixture = fixture;
report.wavenumber = k;
report.quadrature_order = options.QuadratureOrder;
report.tolerance = tol;

switch kind
    case "soft"
        mesh = VolMesh(volPath);
        surface = mesh.boundary();
        g = -exp(1i * k * surface.vtx(:, 3));
        sol = singleLayerDirichletSolve(surface, g, ...
            "Wavenumber", k, "QuadratureOrder", options.QuadratureOrder);
        ref = softSphereScattering(k, 1.0, probes);
        scat = sol.potentialAt(probes);
        report.reference = "soft_sphere_partial_wave_series";
        report.relative_error = ...
            max(abs(scat - ref.scattered) ./ abs(ref.scattered));
        report.series_truncation_tail = ref.truncationTail;
        report.solve_status = string(sol.status);

    case "rigid"
        mesh = VolMesh(volPath);
        surface = mesh.boundary();
        sol = rigidScatteringSolve(surface, "Wavenumber", k, ...
            "QuadratureOrder", options.QuadratureOrder);
        ref = rigidSphereScattering(k, 1.0, probes);
        scat = sol.scatteredAt(probes);
        report.reference = "rigid_sphere_partial_wave_series";
        report.relative_error = ...
            max(abs(scat - ref.scattered) ./ abs(ref.scattered));
        report.condition_number = sol.conditionNumber;
        report.series_truncation_tail = ref.truncationTail;
        report.solve_status = string(sol.status);

    case "coupled_invisible"
        m = FemBemModel(volPath);
        sol = femBemCoupledSolve(m, "Wavenumber", k, ...
            "VolumeSource", 0, "IncidentAmplitude", 1, ...
            "QuadratureOrder", options.QuadratureOrder);
        pinc = exp(1i * k * m.mesh.vtx(:, 3));
        report.reference = "acoustic_invisibility_exact_null";
        report.relative_error = norm(sol.u - pinc) / norm(pinc);
        report.scattered_probe_max = max(abs(sol.exteriorPotentialAt(probes)));
        report.solve_status = string(sol.status);

    case "coupled_fluid"
        k1 = k / 0.7;
        rhor = 1.2;
        m = FemBemModel(volPath);
        sol = femBemCoupledSolve(m, "Wavenumber", k, ...
            "InteriorWavenumber", k1, "DensityRatio", rhor, ...
            "VolumeSource", 0, "IncidentAmplitude", 1, ...
            "QuadratureOrder", options.QuadratureOrder);
        ref = fluidSphereScattering(k, 1.0, m.mesh.vtx, ...
            "InteriorWavenumber", k1, "DensityRatio", rhor);
        report.reference = "anderson_fluid_sphere_series";
        report.interior_wavenumber = k1;
        report.density_ratio = rhor;
        report.relative_error = norm(sol.u - ref.total) / norm(ref.total);
        report.series_truncation_tail = ref.truncationTail;
        report.solve_status = string(sol.status);

    case "focus_adjoint"
        mesh = VolMesh(volPath);
        surface = mesh.boundary();
        nSrc = 8;
        ang = (0:nSrc-1).' / nSrc * 2*pi;
        sources = [2.5*cos(ang), 2.5*sin(ang), -3.0*ones(nSrc, 1)];
        target = [0 0 2.5];
        rng(7);
        amps = randn(nSrc, 1) + 1i*randn(nSrc, 1);
        res = acousticFocusAdjoint(surface, sources, target, k, amps, ...
            "QuadratureOrder", options.QuadratureOrder, "GradientCheck", true);
        report.reference = "adjoint_gradient_vs_finite_difference";
        report.num_sources = nSrc;
        report.relative_error = res.gradientCheckRelError;
        report.forward_affine_residual = res.forwardLinearityResidual;
        report.adjoint_solves = res.adjointSolves;
        report.solve_status = string(res.status);

    case "radiation_force"
        mesh = VolMesh(volPath);
        surface = mesh.boundary();
        sol = rigidScatteringSolve(surface, "Wavenumber", k, ...
            "QuadratureOrder", options.QuadratureOrder);
        t = sol.trace;
        pField = @(X) exp(1i*k*X(:,3)) ...
            + doubleLayerPotentialMatrix(surface, X, k, options.QuadratureOrder) * t;
        rfBem = acousticRadiationForce(pField, k, "ControlRadius", 1.5);
        rfSer = acousticRadiationForce(@(X) rigidSeriesFieldLocal(k, X), k, ...
            "ControlRadius", 1.5);
        report.reference = "brillouin_radiation_stress_bem_vs_series";
        report.force_function_bem = rfBem.forceFunction;
        report.force_function_series = rfSer.forceFunction;
        report.control_radius_residual = rfBem.controlRadiusResidual;
        report.relative_error = abs(rfBem.forceFunction - rfSer.forceFunction) ...
            / abs(rfSer.forceFunction);
        report.solve_status = string(rfBem.status);

    case "fsi"
        model = FemBemModel(volPath);
        probes = [2 0 0; 0 0 3; -1.2 1.6 0];
        % stiff-limit formulation gate: the coupled solve must reproduce
        % the rigid sphere (the interface + BEM coupling is exact,
        % independent of the interior elastic resolution).
        sol = fsiCoupledSolve(model, "Wavenumber", k, ...
            "LongitudinalSpeed", 50, "ShearSpeed", 30, "DensityRatio", 100);
        rigid = rigidSphereScattering(k, 1.0, probes);
        report.reference = "fsi_stiff_limit_vs_rigid_sphere";
        report.relative_error = ...
            max(abs(sol.totalAt(probes) - rigid.total) ./ abs(rigid.total));
        report.solve_status = string(sol.status);

    case "fsi_dtn"
        model = FemBemModel(volPath);
        probes = [2 0 0; 0 0 3; -1.2 1.6 0];
        % same stiff-limit gate on the FAST exterior: the exact spherical
        % Helmholtz DtN / radiating-impedance operator instead of the dense
        % Galerkin BEM - reproduces the rigid sphere with no dense N^2
        % assembly (fail-loud if the truncation is not a sphere). Acoustic
        % waves do not use the Kelvin-boundary label in this lab policy.
        sol = fsiCoupledSolve(model, "Wavenumber", k, ...
            "LongitudinalSpeed", 50, "ShearSpeed", 30, "DensityRatio", 100, ...
            "ExteriorMethod", "dtn");
        rigid = rigidSphereScattering(k, 1.0, probes);
        report.reference = "fsi_dtn_stiff_limit_vs_rigid_sphere";
        report.exterior_method = string(sol.exteriorMethod);
        report.dtn_degree = sol.dtn.degree;
        report.dtn_num_modes = sol.dtn.numModes;
        report.relative_error = ...
            max(abs(sol.totalAt(probes) - rigid.total) ./ abs(rigid.total));
        report.solve_status = string(sol.status);

    case "thrust_adjoint"
        % the wavefront-synthesis thrust adjoint through the FSI solve: the
        % radiation force on an ELASTIC bead and its gradient wrt the phased-
        % array amplitudes. relative_error = max(form-vs-direct, gradient-vs-FD).
        model = FemBemModel(volPath);
        nSrc = 4;
        ang = (0:nSrc-1).' / nSrc * 2*pi;
        srcs = [2.5*cos(ang), 2.5*sin(ang), -3*ones(nSrc, 1)];   % ring below the bead
        amps = exp(1i * ang * 0.7);
        adj = elasticThrustAdjoint(model, srcs, k, amps, ...
            "NMu", 8, "NPhi", 12, "GradientCheck", true);
        report.reference = "thrust_adjoint_form_vs_direct_and_gradient_vs_fd";
        report.force_z = adj.force(3);
        report.consistency_error = adj.consistencyError;
        report.gradient_check = adj.gradientCheckRelError;
        report.control_radius_residual = adj.controlRadiusResidual;
        report.relative_error = max(adj.consistencyError, adj.gradientCheckRelError);
        report.solve_status = string(adj.status);
end

report.pass = isfinite(report.relative_error) ...
    && report.relative_error <= tol ...
    && report.solve_status == "ok";
if report.pass
    report.status = "ok";
else
    report.status = "needs_attention";
end
end


function p = rigidSeriesFieldLocal(k, X)
%RIGIDSERIESFIELDLOCAL Total field of a rigid unit sphere, +z plane wave.
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
        dj = sphjLocal(l-1, k) - (l+1)/k * sphjLocal(l, k);
        dh = sphhLocal(l-1, k) - (l+1)/k * sphhLocal(l, k);
        sc = sc - (1i^l)*(2*l+1)*dj/dh * sphhLocal(l, k*r) * Pl;
    end
    p(m) = exp(1i*k*x(3)) + sc;
end
end

function j = sphjLocal(l, x)
j = sqrt(pi/(2*x)) * besselj(l+0.5, x);
end

function h = sphhLocal(l, x)
h = sqrt(pi/(2*x)) * (besselj(l+0.5, x) + 1i*bessely(l+0.5, x));
end
