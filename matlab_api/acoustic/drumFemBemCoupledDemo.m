function scene = drumFemBemCoupledDemo(options)
%DRUMFEMBEMCOUPLEDDEMO Reduced FEM/BEM drum acoustics teaching demo.
%
%   scene = drumFemBemCoupledDemo() builds a readable time-domain
%   vibro-acoustic demo for a real drum-like object:
%
%       top membrane FEM mode  -> internal cavity pressure FEM mode
%       bottom membrane FEM mode + shell leakage
%       -> exterior retarded-potential boundary integrals (BEM view)
%
%   This is deliberately a reduced-order teaching model, not a production
%   drum solver.  Its purpose is to make the coupling topology visible before
%   moving to P1 volume FEM + P1 BEM on a .vol mesh.

arguments
    options.Radius (1,1) double {mustBePositive} = 0.10
    options.Depth (1,1) double {mustBePositive} = 0.075
    options.RimThickness (1,1) double {mustBeNonnegative} = 0.012
    options.BoundaryRadius (1,1) double {mustBePositive} = 0.28
    options.SoundSpeed (1,1) double {mustBePositive} = 343.0
    options.Density (1,1) double {mustBePositive} = 1.2
    options.TopFrequency (1,1) double {mustBePositive} = 220.0
    options.BottomFrequency (1,1) double {mustBePositive} = 185.0
    options.CavityFrequency (1,1) double {mustBePositive} = 620.0
    options.MembraneDamping (1,1) double {mustBeNonnegative} = 0.035
    options.CavityDamping (1,1) double {mustBeNonnegative} = 0.055
    options.CouplingStrength (1,1) double {mustBeNonnegative} = 1.0e4
    options.CavityGain (1,1) double {mustBeNonnegative} = 5.0e5
    options.SideLeakage (1,1) double {mustBeNonnegative} = 2.5e3
    options.CavityDisplayGain (1,1) double {mustBeNonnegative} = 0.0
    options.ForcePerModalMass (1,1) double = 1.0
    options.Excitation (1,1) string {mustBeMember(options.Excitation, ["impact", "step"])} = "impact"
    options.ImpactDuration (1,1) double {mustBePositive} = 1.8e-4
    options.NumX (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumX, 24)} = 120
    options.NumZ (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumZ, 24)} = 120
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 8)} = 72
    options.TMax (1,1) double {mustBePositive} = 2.4e-3
    options.NumSourceRadial (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumSourceRadial, 2)} = 10
    options.NumSourceAzimuth (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumSourceAzimuth, 3)} = 20
    options.NumSideAxial (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumSideAxial, 2)} = 8
end

a = options.Radius;
depth = options.Depth;
boundaryRadius = options.BoundaryRadius;
if boundaryRadius <= sqrt(a^2 + depth^2)
    error("drumFemBemCoupledDemo:boundaryRadius", ...
        "BoundaryRadius must enclose the drum body.");
end

t = linspace(0, options.TMax, options.NumTime);
motion = solveReducedFemDrum(t, options);

x = linspace(-boundaryRadius, boundaryRadius, options.NumX);
z = linspace(-boundaryRadius, boundaryRadius, options.NumZ);
[X, Z] = meshgrid(x, z);
rho = hypot(X, Z);
domain = rho <= boundaryRadius;

[masks, geometry] = drumMasks(X, Z, options);
pressure = zeros(numel(z), numel(x), numel(t));

exteriorObservation = domain & ~masks.cavity & ~masks.drum_frame;
upperExterior = exteriorObservation & Z > 0;
lowerExterior = exteriorObservation & Z < -depth;
lateralExterior = exteriorObservation & Z < 0 & Z > -depth;

upperVec = upperExterior(exteriorObservation);
lowerVec = lowerExterior(exteriorObservation);
lateralVec = lateralExterior(exteriorObservation);
componentMax = struct( ...
    "top_to_lateral_exterior", 0, ...
    "top_to_lower_exterior", 0, ...
    "bottom_to_upper_exterior", 0, ...
    "side_to_upper_exterior", 0, ...
    "side_to_lower_exterior", 0);
for k = 1:numel(t)
    p = zeros(size(X));
    pTop = diskRayleighAt( ...
        X(exteriorObservation), Z(exteriorObservation), t(k), t, motion.top_acceleration, ...
        0.0, +1, options);
    pBottom = diskRayleighAt( ...
        X(exteriorObservation), Z(exteriorObservation), t(k), t, motion.bottom_acceleration, ...
        -depth, -1, options);
    pSide = sideLeakageAt( ...
        X(exteriorObservation), Z(exteriorObservation), t(k), t, motion.side_leakage_acceleration, ...
        options);
    p(exteriorObservation) = pTop + pBottom + pSide;

    componentMax.top_to_lateral_exterior = maxMaskedComponent( ...
        componentMax.top_to_lateral_exterior, pTop, lateralVec);
    componentMax.top_to_lower_exterior = maxMaskedComponent( ...
        componentMax.top_to_lower_exterior, pTop, lowerVec);
    componentMax.bottom_to_upper_exterior = maxMaskedComponent( ...
        componentMax.bottom_to_upper_exterior, pBottom, upperVec);
    componentMax.side_to_upper_exterior = maxMaskedComponent( ...
        componentMax.side_to_upper_exterior, pSide, upperVec);
    componentMax.side_to_lower_exterior = maxMaskedComponent( ...
        componentMax.side_to_lower_exterior, pSide, lowerVec);

    if options.CavityDisplayGain > 0
        p(masks.cavity) = options.CavityDisplayGain * motion.cavity_pressure(k);
    end
    pressure(:, :, k) = p;
end

scene = struct();
scene.kind = "drum_reduced_fem_bem_coupled_scene";
scene.x = x;
scene.z = z;
scene.t = t;
scene.time_indices = 1:numel(t);
scene.pressure = pressure;
scene.motion = motion;
scene.source = "reduced FEM membrane/cavity model coupled to exterior BEM-style retarded integrals";
scene.boundary_type = "full spherical high-order impedance absorbing boundary visualization";
scene.coupling = struct();
scene.coupling.kind = "reduced_fem_internal_air_plus_exterior_bem";
scene.coupling.fem_dofs = ["top_membrane_mode", "bottom_membrane_mode", "internal_cavity_pressure"];
scene.coupling.bem_surfaces = ["top_head_exterior", "bottom_head_exterior", "side_shell_leakage"];
scene.coupling.note = "Every BEM boundary source is evaluated at every exterior observation point with the same retarded free-space Green kernel. The internal cavity pressure drives the coupling but is not color-mapped by default.";
scene.coupling.cavity_display_gain = options.CavityDisplayGain;
scene.bem = struct();
scene.bem.kernel = "causal retarded free-space Green function";
scene.bem.observation_rule = "top, bottom, and side boundary sources all contribute to every exterior air observation point";
scene.bem.unknown_view = "FEM modal accelerations provide Neumann boundary data for the BEM teaching layer";
scene.visualization = struct();
scene.visualization.field = "propagating_air_pressure_only";
scene.visualization.internal_cavity = "special coupling state, not a plotted pressure field";
scene.geometry = geometry;
scene.masks = masks;
scene.axis = struct();
scene.axis.equal = true;
scene.axis.x_limits = [x(1), x(end)];
scene.axis.z_limits = [z(1), z(end)];
scene.axis.pixel_size = [numel(x), numel(z)];
scene.summary = struct();
scene.summary.max_abs_pressure = max(abs(pressure), [], "all");
scene.summary.top_peak_acceleration = max(abs(motion.top_acceleration));
scene.summary.bottom_peak_acceleration = max(abs(motion.bottom_acceleration));
scene.summary.cavity_peak_pressure = max(abs(motion.cavity_pressure));
lowerMask3 = repmat(Z < -depth, 1, 1, numel(t));
scene.summary.lower_half_wave_present = max(abs(pressure(lowerMask3)), [], "all") > 0;
scene.summary.internal_resonance_present = scene.summary.cavity_peak_pressure > 0;
scene.summary.bem_cross_direction = componentMax;
cavityMask3 = repmat(masks.cavity, 1, 1, numel(t));
scene.summary.cavity_display_peak = max(abs(pressure(cavityMask3)), [], "all");
scene.checks = struct();
scene.checks.finite_pressure = all(isfinite(pressure), "all");
scene.checks.lower_half_wave_present = scene.summary.lower_half_wave_present;
scene.checks.internal_cavity_coupled = scene.summary.internal_resonance_present;
scene.checks.cavity_not_colormapped = scene.summary.cavity_display_peak == 0;
crossScale = max(scene.summary.max_abs_pressure, eps);
scene.checks.top_source_reaches_lateral_exterior = ...
    componentMax.top_to_lateral_exterior > 1e-10 * crossScale;
scene.checks.top_source_reaches_lower_exterior = ...
    componentMax.top_to_lower_exterior > 1e-10 * crossScale;
scene.checks.bottom_source_reaches_upper_exterior = ...
    componentMax.bottom_to_upper_exterior > 1e-10 * crossScale;
scene.checks.side_source_reaches_upper_exterior = ...
    componentMax.side_to_upper_exterior > 1e-10 * crossScale;
scene.checks.side_source_reaches_lower_exterior = ...
    componentMax.side_to_lower_exterior > 1e-10 * crossScale;
scene.checks.axis_equal_grid = numel(x) == numel(z) && ...
    abs((x(end) - x(1)) - (z(end) - z(1))) < 1e-12;
if all(structfun(@(v) logical(v), scene.checks))
    scene.status = "ok";
else
    scene.status = "needs_attention";
end
end


function current = maxMaskedComponent(current, values, mask)
if any(mask)
    current = max(current, max(abs(values(mask)), [], "all"));
end
end


function motion = solveReducedFemDrum(t, options)
wTop = 2 * pi * options.TopFrequency;
wBottom = 2 * pi * options.BottomFrequency;
wCavity = 2 * pi * options.CavityFrequency;
zetaM = min(options.MembraneDamping, 0.999);
zetaC = min(options.CavityDamping, 0.999);
alpha = options.CouplingStrength;
beta = options.CavityGain;
force = options.ForcePerModalMass;
forceAt = @(time) impactForce(time, force, options.Excitation, options.ImpactDuration);

opts = odeset("RelTol", 1e-7, "AbsTol", 1e-9);
odeRhs = @(time, y) [
    y(2)
    forceAt(time) - 2*zetaM*wTop*y(2) - wTop^2*y(1) - alpha*y(5)
    y(4)
    alpha*y(5) - 2*zetaM*wBottom*y(4) - wBottom^2*y(3)
    y(6)
    beta*(y(2) - y(4)) - 2*zetaC*wCavity*y(6) - wCavity^2*y(5)
    ];
[~, y] = ode45(odeRhs, t, zeros(6, 1), opts);

forceHistory = forceAt(t(:));
topAcceleration = forceHistory - 2*zetaM*wTop*y(:, 2) - wTop^2*y(:, 1) - alpha*y(:, 5);
bottomAcceleration = alpha*y(:, 5) - 2*zetaM*wBottom*y(:, 4) - wBottom^2*y(:, 3);
sideLeakageAcceleration = options.SideLeakage * y(:, 5);

motion = struct();
motion.kind = "reduced_fem_membrane_cavity_motion";
motion.t = t;
motion.top_displacement = y(:, 1).';
motion.top_velocity = y(:, 2).';
motion.bottom_displacement = y(:, 3).';
motion.bottom_velocity = y(:, 4).';
motion.cavity_pressure = y(:, 5).';
motion.cavity_pressure_rate = y(:, 6).';
motion.force = forceHistory.';
motion.top_acceleration = topAcceleration.';
motion.bottom_acceleration = bottomAcceleration.';
motion.side_leakage_acceleration = sideLeakageAcceleration.';
motion.parameters = struct( ...
    "top_frequency", options.TopFrequency, ...
    "bottom_frequency", options.BottomFrequency, ...
    "cavity_frequency", options.CavityFrequency, ...
    "coupling_strength", options.CouplingStrength, ...
    "cavity_gain", options.CavityGain, ...
    "side_leakage", options.SideLeakage, ...
    "excitation", options.Excitation, ...
    "impact_duration", options.ImpactDuration);
end


function f = impactForce(t, amplitude, excitation, duration)
switch excitation
    case "step"
        f = amplitude * double(t >= 0);
    otherwise
        f = zeros(size(t));
        active = t >= 0 & t <= duration;
        f(active) = amplitude * sin(pi * t(active) / duration);
end
end


function value = diskRayleighAt(x, z, tk, t, acceleration, zSource, normalSign, options)
[srcR, srcTheta, srcWeight] = diskPolarQuadrature(options.Radius, ...
    options.NumSourceRadial, options.NumSourceAzimuth);
mode = drumMode01(srcR / options.Radius);
sourceWeight = srcWeight .* mode;
value = zeros(size(x));
for m = 1:numel(srcR)
    xs = srcR(m) * cos(srcTheta(m));
    ys = srcR(m) * sin(srcTheta(m));
    dist = sqrt((x - xs).^2 + ys.^2 + (z - zSource).^2);
    tau = tk - dist / options.SoundSpeed;
    delayed = interp1(t, acceleration, tau, "linear", 0);
    value = value + normalSign * sourceWeight(m) * delayed ./ max(dist, eps);
end
value = options.Density / (2 * pi) * value;
end


function value = sideLeakageAt(x, z, tk, t, acceleration, options)
radius = options.Radius + options.RimThickness;
zSamples = linspace(-options.Depth, 0, options.NumSideAxial);
theta = linspace(0, 2*pi, options.NumSourceAzimuth + 1);
theta(end) = [];
dz = options.Depth / max(options.NumSideAxial - 1, 1);
dtheta = 2 * pi / options.NumSourceAzimuth;
weight = radius * dz * dtheta;
value = zeros(size(x));
for iz = 1:numel(zSamples)
    for it = 1:numel(theta)
        xs = radius * cos(theta(it));
        ys = radius * sin(theta(it));
        zs = zSamples(iz);
        dist = sqrt((x - xs).^2 + ys.^2 + (z - zs).^2);
        tau = tk - dist / options.SoundSpeed;
        delayed = interp1(t, acceleration, tau, "linear", 0);
        value = value + weight * delayed ./ max(dist, eps);
    end
end
value = options.Density / (4 * pi) * value;
end


function [masks, geometry] = drumMasks(X, Z, options)
a = options.Radius;
depth = options.Depth;
outerRadius = a + options.RimThickness;
boundaryRadius = options.BoundaryRadius;
dx = abs(X(1, 2) - X(1, 1));
dz = abs(Z(2, 1) - Z(1, 1));
tol = 0.75 * max(dx, dz);

cavity = abs(X) < a & Z < 0 & Z > -depth;
insideCylinder = abs(X) <= outerRadius & Z >= -depth & Z <= 0;
hollow = abs(X) < a & Z > -depth & Z < 0;
drumFrame = insideCylinder & ~hollow;

sideWall = abs(abs(X) - outerRadius) <= tol & Z >= -depth & Z <= 0;
topRim = abs(Z) <= tol & abs(X) >= a & abs(X) <= outerRadius;
bottomRim = abs(Z + depth) <= tol & abs(X) >= a & abs(X) <= outerRadius;
frameOutline = sideWall | topRim | bottomRim;
topMembrane = abs(Z) <= tol & abs(X) <= a;
bottomMembrane = abs(Z + depth) <= tol & abs(X) <= a;

rho = hypot(X, Z);
boundary = abs(rho - boundaryRadius) <= tol;
boundaryDomain = rho <= boundaryRadius;

masks = struct();
masks.boundary_domain = boundaryDomain;
masks.high_order_impedance_boundary = boundary;
masks.drum_frame = drumFrame;
masks.frame_outline = frameOutline;
masks.membrane = topMembrane | bottomMembrane;
masks.top_membrane = topMembrane;
masks.bottom_membrane = bottomMembrane;
masks.cavity = cavity;

geometry = struct();
geometry.radius = a;
geometry.outer_radius = outerRadius;
geometry.depth = depth;
geometry.frame_depth = depth;
geometry.boundary_radius = boundaryRadius;
geometry.struck_surface = "top membrane at z=0";
geometry.radiation_model = "two-head drum with internal cavity and side leakage";
geometry.description = "cylindrical drum body inside a full spherical high-order impedance boundary";
end


function [r, theta, w] = diskPolarQuadrature(radius, nr, nth)
dr = radius / nr;
rCenters = ((1:nr).' - 0.5) * dr;
dtheta = 2 * pi / nth;
thetaCenters = ((1:nth) - 0.5) * dtheta;
[R, T] = ndgrid(rCenters, thetaCenters);
r = R(:);
theta = T(:);
w = r * dr * dtheta;
end


function y = drumMode01(rho)
alpha01 = 2.4048255577;
y = besselj(0, alpha01 * rho);
end
