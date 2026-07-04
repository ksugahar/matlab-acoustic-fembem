function info = writeDrumHighOrderImpedanceGif(scene, gifPath, options)
%WRITEDRUMHIGHORDERIMPEDANCEGIF Write a drum wave GIF with boundary overlays.
%
%   info = writeDrumHighOrderImpedanceGif(scene, gifPath) writes an animated
%   GIF from drumHighOrderImpedanceScene.  The pressure map is direct
%   indexed-image output for batch/headless runs; drum frame, membrane, and
%   high-order impedance boundary are overlaid as image masks.

arguments
    scene (1,1) struct
    gifPath (1,1) string
    options.DelayTime (1,1) double {mustBePositive} = 0.04
    options.LoopCount (1,1) double {mustBeNonnegative} = Inf
    options.Colormap (1,1) string = "turbo"
    options.NumFieldColors (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumFieldColors, 16), mustBeLessThanOrEqual(options.NumFieldColors, 248)} = 240
    options.FlipVertical (1,1) logical = true
end

validateScene(scene);
validateLoopCount(options.LoopCount);

gifPath = string(gifPath);
outDir = fileparts(gifPath);
if strlength(outDir) > 0 && ~isfolder(outDir)
    mkdir(outDir);
end

[map, idx] = palette(options.Colormap, options.NumFieldColors);
scale = max(abs(scene.pressure), [], "all");
if ~isfinite(scale) || scale <= 0
    scale = 1;
end

for k = 1:size(scene.pressure, 3)
    frame = scene.pressure(:, :, k);
    image = pressureToIndexed(frame, scale, options.NumFieldColors);
    image(~scene.masks.boundary_domain) = idx.background;
    image(scene.masks.drum_frame) = idx.frame_fill;
    image(scene.masks.high_order_impedance_boundary) = idx.boundary;
    image(scene.masks.frame_outline) = idx.frame_outline;
    image(scene.masks.membrane) = idx.membrane;
    if options.FlipVertical
        image = flipud(image);
    end

    if k == 1
        imwrite(image, map, gifPath, "gif", ...
            "LoopCount", options.LoopCount, ...
            "DelayTime", options.DelayTime);
    else
        imwrite(image, map, gifPath, "gif", ...
            "WriteMode", "append", ...
            "DelayTime", options.DelayTime);
    end
end

info = struct();
info.kind = "drum_high_order_impedance_gif";
info.path = gifPath;
info.num_frames = size(scene.pressure, 3);
info.delay_time = options.DelayTime;
info.loop_count = options.LoopCount;
info.pressure_scale = scale;
info.boundary_type = scene.boundary_type;
info.boundary_radius = scene.geometry.boundary_radius;
info.frame_depth = scene.geometry.frame_depth;
if isfield(scene, "axis")
    info.axis_equal = scene.axis.equal;
end
info.flip_vertical = options.FlipVertical;
end


function validateScene(scene)
required = ["pressure", "masks", "geometry", "boundary_type"];
for name = required
    if ~isfield(scene, name)
        error("writeDrumHighOrderImpedanceGif:InvalidScene", ...
            "scene.%s is required.", name);
    end
end
maskNames = ["boundary_domain", "high_order_impedance_boundary", ...
    "drum_frame", "frame_outline", "membrane"];
for name = maskNames
    if ~isfield(scene.masks, name)
        error("writeDrumHighOrderImpedanceGif:InvalidScene", ...
            "scene.masks.%s is required.", name);
    end
end
end


function validateLoopCount(loopCount)
if isfinite(loopCount) && loopCount ~= floor(loopCount)
    error("writeDrumHighOrderImpedanceGif:InvalidLoopCount", ...
        "LoopCount must be an integer or Inf.");
end
end


function [map, idx] = palette(name, n)
try
    fieldMap = feval(name, n);
catch
    fieldMap = parula(n);
end
map = [
    fieldMap
    0.96 0.96 0.94   % background
    0.33 0.33 0.34   % drum frame fill
    0.05 0.05 0.05   % frame outline
    1.00 0.85 0.16   % membrane
    0.00 0.00 0.00   % high-order impedance boundary
    ];

idx = struct();
idx.background = uint8(n);
idx.frame_fill = uint8(n + 1);
idx.frame_outline = uint8(n + 2);
idx.membrane = uint8(n + 3);
idx.boundary = uint8(n + 4);
end


function image = pressureToIndexed(frame, scale, numColors)
value = max(-1, min(1, frame / scale));
image = uint8(round((value + 1) * 0.5 * (numColors - 1)));
end
