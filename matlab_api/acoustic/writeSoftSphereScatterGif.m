function info = writeSoftSphereScatterGif(field, gifPath, options)
%WRITESOFTSPHERESCATTERGIF Write an animated GIF of soft-sphere pulse scattering.
%
%   info = writeSoftSphereScatterGif(field, gifPath) writes the x-z total-field
%   snapshots from softSphereScatterField to an animated GIF.  The pressure array
%   is mapped directly to an indexed image (points inside the sphere and outside
%   the colour clip are drawn in a fixed sphere/background colour), so the writer
%   works in headless MATLAB runs without a visible figure window.

arguments
    field (1,1) struct
    gifPath (1,1) string
    options.TimeIndices = []
    options.DelayTime (1,1) double {mustBePositive} = 0.18
    options.LoopCount (1,1) double {mustBeNonnegative} = Inf
    options.ColorScale double = []
    options.Colormap (1,1) string = "turbo"
    options.NumColors (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumColors,1), mustBeLessThanOrEqual(options.NumColors,256)} = 254
end

validateField(field);

if isempty(options.TimeIndices)
    timeIndices = 1:numel(field.time);
else
    timeIndices = double(options.TimeIndices(:)).';
    mustBeInteger(timeIndices);
    mustBePositive(timeIndices);
    timeIndices = min(timeIndices, numel(field.time));
end

gifPath = string(gifPath);
outDir = fileparts(gifPath);
if strlength(outDir) > 0 && ~isfolder(outDir)
    mkdir(outDir);
end

cmap = makeColormap(options.Colormap, options.NumColors);
sphereColor = [0.15 0.15 0.15];                  % index NumColors+1
cmapFull = [cmap; sphereColor];

scale = options.ColorScale;
if isempty(scale)
    scale = prctile(abs(field.pressure(~isnan(field.pressure))), 97);
end
if ~isfinite(scale) || scale <= 0
    scale = 1;
end

for k = 1:numel(timeIndices)
    frame = field.pressure(:, :, timeIndices(k));         % [nz, nx] (imagesc layout)
    indexed = pressureToIndexed(frame, scale, options.NumColors);
    if k == 1
        imwrite(indexed, cmapFull, gifPath, "gif", ...
            "LoopCount", options.LoopCount, "DelayTime", options.DelayTime);
    else
        imwrite(indexed, cmapFull, gifPath, "gif", ...
            "WriteMode", "append", "DelayTime", options.DelayTime);
    end
end

info = struct();
info.kind = "soft_sphere_scatter_gif";
info.path = gifPath;
info.num_frames = numel(timeIndices);
info.time_indices = timeIndices;
info.delay_time = options.DelayTime;
info.loop_count = options.LoopCount;
info.colormap = options.Colormap;
info.num_colors = options.NumColors;
info.color_scale = scale;
info.time_range = [field.time(timeIndices(1)), field.time(timeIndices(end))];
end


function validateField(field)
required = ["pressure", "x", "z", "time"];
for name = required
    if ~isfield(field, name)
        error("writeSoftSphereScatterGif:InvalidField", "field.%s is required.", name);
    end
end
if ndims(field.pressure) ~= 3
    error("writeSoftSphereScatterGif:InvalidField", ...
        "field.pressure must be a 3-D array [nz, nx, nt].");
end
if size(field.pressure, 3) ~= numel(field.time)
    error("writeSoftSphereScatterGif:InvalidField", ...
        "field.pressure third dimension must match numel(field.time).");
end
end


function cmap = makeColormap(name, n)
try
    cmap = feval(name, n);
catch
    cmap = parula(n);
end
end


function indexed = pressureToIndexed(frame, scale, numColors)
% NaN (inside sphere / masked) -> the extra sphere-colour index (numColors).
value = max(-1, min(1, frame / scale));
idx = uint8(round((value + 1) * 0.5 * (numColors - 1)));
idx(isnan(frame)) = numColors;                    % 0-based extra colour row
indexed = idx;
end
