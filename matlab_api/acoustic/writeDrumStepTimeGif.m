function info = writeDrumStepTimeGif(field, gifPath, options)
%WRITEDRUMSTEPTIMEGIF Write an animated GIF of a drum pressure time field.
%
%   info = writeDrumStepTimeGif(field, gifPath) writes the r-z pressure
%   snapshots from drumStepTimeField to an animated GIF.  The writer maps
%   the pressure array directly to an indexed image, so it works in
%   headless MATLAB runs and does not depend on a visible figure window.

arguments
    field (1,1) struct
    gifPath (1,1) string
    options.TimeIndices = []
    options.DelayTime (1,1) double {mustBePositive} = 0.04
    options.LoopCount (1,1) double {mustBeNonnegative} = Inf
    options.Normalize (1,1) logical = true
    options.Colormap (1,1) string = "turbo"
    options.NumColors (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumColors, 1), mustBeLessThanOrEqual(options.NumColors, 256)} = 256
end

validateField(field);
validateLoopCount(options.LoopCount);

if isempty(options.TimeIndices)
    timeIndices = 1:numel(field.t);
else
    timeIndices = double(options.TimeIndices(:)).';
    mustBeInteger(timeIndices);
    mustBePositive(timeIndices);
    timeIndices = min(timeIndices, numel(field.t));
end

gifPath = string(gifPath);
outDir = fileparts(gifPath);
if strlength(outDir) > 0 && ~isfolder(outDir)
    mkdir(outDir);
end

cmap = makeColormap(options.Colormap, options.NumColors);
scale = max(abs(field.pressure), [], "all");
if ~isfinite(scale) || scale <= 0
    scale = 1;
end

for k = 1:numel(timeIndices)
    frame = field.pressure(:, :, timeIndices(k)).';
    indexed = pressureToIndexed(frame, scale, options.Normalize, options.NumColors);
    if k == 1
        imwrite(indexed, cmap, gifPath, "gif", ...
            "LoopCount", options.LoopCount, ...
            "DelayTime", options.DelayTime);
    else
        imwrite(indexed, cmap, gifPath, "gif", ...
            "WriteMode", "append", ...
            "DelayTime", options.DelayTime);
    end
end

info = struct();
info.kind = "drum_step_time_field_gif";
info.path = gifPath;
info.num_frames = numel(timeIndices);
info.time_indices = timeIndices;
info.delay_time = options.DelayTime;
info.loop_count = options.LoopCount;
info.normalize = options.Normalize;
info.colormap = options.Colormap;
info.num_colors = options.NumColors;
info.pressure_scale = scale;
info.time_range = [field.t(timeIndices(1)), field.t(timeIndices(end))];
end


function validateLoopCount(loopCount)
if isfinite(loopCount) && loopCount ~= floor(loopCount)
    error("writeDrumStepTimeGif:InvalidLoopCount", ...
        "LoopCount must be an integer or Inf.");
end
end


function validateField(field)
required = ["pressure", "r", "z", "t"];
for name = required
    if ~isfield(field, name)
        error("writeDrumStepTimeGif:InvalidField", ...
            "field.%s is required.", name);
    end
end
if ndims(field.pressure) ~= 3
    error("writeDrumStepTimeGif:InvalidField", ...
        "field.pressure must be a 3-D array [nr, nz, nt].");
end
if size(field.pressure, 3) ~= numel(field.t)
    error("writeDrumStepTimeGif:InvalidField", ...
        "field.pressure third dimension must match numel(field.t).");
end
end


function cmap = makeColormap(name, n)
try
    cmap = feval(name, n);
catch
    cmap = parula(n);
end
end


function indexed = pressureToIndexed(frame, scale, normalize, numColors)
if normalize
    value = frame / scale;
else
    value = frame;
    maxAbs = max(abs(value), [], "all");
    if isfinite(maxAbs) && maxAbs > 0
        value = value / maxAbs;
    end
end
value = max(-1, min(1, value));
indexed = uint8(round((value + 1) * 0.5 * (numColors - 1)));
end
