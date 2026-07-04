function scene = drumHighOrderImpedanceScene(field, options)
%DRUMHIGHORDERIMPEDANCESCENE Build a drum + high-order impedance boundary scene.
%
%   scene = drumHighOrderImpedanceScene(field) converts the axisymmetric
%   top-side pressure returned by drumStepTimeField into an x-z cut-plane
%   scene.  The scene contains a cylindrical drum body, its struck top
%   membrane, and a full spherical truncation boundary annotated as the
%   high-order impedance / absorbing-boundary lane used for Radia-style open
%   boundary experiments.  This is a visualization scaffold; the pressure is
%   still the readable Rayleigh time-domain model from drumStepTimeField.

arguments
    field (1,1) struct
    options.BoundaryRadius (1,1) double {mustBeNonnegative} = 0
    options.FrameDepth (1,1) double {mustBeNonnegative} = 0
    options.RimThickness (1,1) double {mustBeNonnegative} = 0
    options.NumX (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumX, 8)} = 150
    options.NumZ (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumZ, 8)} = 150
    options.AxisEqual (1,1) logical = true
    options.TimeIndices = []
end

validateField(field);

radius = field.radius;
boundaryRadius = options.BoundaryRadius;
if boundaryRadius <= 0
    boundaryRadius = 0.94 * min(field.r(end), field.z(end));
end
boundaryRadius = min(boundaryRadius, min(field.r(end), field.z(end)));

frameDepth = options.FrameDepth;
if frameDepth <= 0
    frameDepth = 0.45 * radius;
end

rimThickness = options.RimThickness;
if rimThickness <= 0
    rimThickness = 0.13 * radius;
end
outerRadius = radius + rimThickness;

numX = options.NumX;
numZ = options.NumZ;
if options.AxisEqual
    numZ = numX;
end

if isempty(options.TimeIndices)
    timeIndices = 1:numel(field.t);
else
    timeIndices = double(options.TimeIndices(:)).';
    mustBeInteger(timeIndices);
    mustBePositive(timeIndices);
    timeIndices = min(timeIndices, numel(field.t));
end

x = linspace(-boundaryRadius, boundaryRadius, numX);
z = linspace(-boundaryRadius, boundaryRadius, numZ);
[X, Z] = meshgrid(x, z);
R = abs(X);
rho = hypot(X, Z);

domain = rho <= boundaryRadius;
topRadiatingDomain = domain & Z >= field.z(1);
interpolationDomain = topRadiatingDomain & R <= field.r(end) & Z <= field.z(end);
pressure = zeros(numel(z), numel(x), numel(timeIndices));
for k = 1:numel(timeIndices)
    pressure(:, :, k) = interpolatePressure(field, timeIndices(k), ...
        R, Z, interpolationDomain);
end

[masks, geometry] = buildMasks(X, Z, radius, outerRadius, ...
    frameDepth, boundaryRadius);

scene = struct();
scene.kind = "drum_high_order_impedance_scene";
scene.x = x;
scene.z = z;
scene.t = field.t(timeIndices);
scene.time_indices = timeIndices;
scene.pressure = pressure;
scene.source = "drumStepTimeField Rayleigh pressure from the struck top membrane";
scene.boundary_type = "full spherical high-order impedance absorbing boundary visualization";
scene.boundary_notes = [
    "The display boundary is a full spherical truncation surface."
    "The drum top surface at z=0 is struck; the cylinder body sits below it."
    "It is intentionally not labelled Kelvin; use it as the high-order"
    "impedance/ABC lane until the Radia production name is finalized."
    ];
scene.geometry = geometry;
scene.masks = masks;
scene.axis = struct();
scene.axis.equal = options.AxisEqual;
scene.axis.x_limits = [x(1), x(end)];
scene.axis.z_limits = [z(1), z(end)];
scene.axis.pixel_size = [numel(x), numel(z)];
scene.summary = struct();
scene.summary.num_frames = numel(timeIndices);
scene.summary.max_abs_pressure = max(abs(pressure), [], "all");
scene.summary.boundary_radius = boundaryRadius;
scene.summary.frame_depth = frameDepth;
scene.summary.outer_radius = outerRadius;
end


function validateField(field)
required = ["pressure", "r", "z", "t", "radius"];
for name = required
    if ~isfield(field, name)
        error("drumHighOrderImpedanceScene:InvalidField", ...
            "field.%s is required.", name);
    end
end
if ndims(field.pressure) ~= 3
    error("drumHighOrderImpedanceScene:InvalidField", ...
        "field.pressure must be a 3-D array [nr, nz, nt].");
end
end


function pressure = interpolatePressure(field, timeIndex, R, Z, mask)
pressure = zeros(size(R));
if ~any(mask, "all")
    return
end
interp = griddedInterpolant({field.r, field.z}, ...
    field.pressure(:, :, timeIndex), "linear", "none");
values = interp(R(mask), Z(mask));
values(~isfinite(values)) = 0;
pressure(mask) = values;
end


function [masks, geometry] = buildMasks(X, Z, radius, outerRadius, ...
        frameDepth, boundaryRadius)
dx = abs(X(1, 2) - X(1, 1));
dz = abs(Z(2, 1) - Z(1, 1));
tol = 0.75 * max(dx, dz);

insideCylinder = abs(X) <= outerRadius & Z >= -frameDepth & Z <= 0;
hollow = abs(X) < radius & Z > -frameDepth & Z < 0;
drumFrame = insideCylinder & ~hollow;

sideWall = abs(abs(X) - outerRadius) <= tol & Z >= -frameDepth & Z <= 0;
backWall = abs(Z + frameDepth) <= tol & abs(X) <= outerRadius;
frontRim = abs(Z) <= tol & abs(X) >= radius & abs(X) <= outerRadius;
frameOutline = sideWall | backWall | frontRim;

membrane = abs(Z) <= tol & abs(X) <= radius;

rho = hypot(X, Z);
boundary = abs(rho - boundaryRadius) <= tol;
boundaryDomain = rho <= boundaryRadius;

masks = struct();
masks.boundary_domain = boundaryDomain;
masks.high_order_impedance_boundary = boundary;
masks.drum_frame = drumFrame;
masks.frame_outline = frameOutline;
masks.membrane = membrane;

geometry = struct();
geometry.radius = radius;
geometry.outer_radius = outerRadius;
geometry.frame_depth = frameDepth;
geometry.boundary_radius = boundaryRadius;
geometry.struck_surface = "top membrane at z=0";
geometry.description = "cylindrical drum cross-section inside a full spherical high-order impedance boundary";
end
