function ax = plotDrumStepTimeField(field, timeIndex, options)
%PLOTDRUMSTEPTIMEFIELD Plot an r-z pressure snapshot from drumStepTimeField.
%
%   ax = plotDrumStepTimeField(field, timeIndex) draws a normalized pressure
%   map.  The horizontal axis is radius from the drum center; the vertical
%   axis is distance away from the drumhead.

arguments
    field (1,1) struct
    timeIndex (1,1) double {mustBeInteger, mustBePositive} = 1
    options.Parent = []
    options.Normalize (1,1) logical = true
end

timeIndex = min(timeIndex, numel(field.t));
snapshot = field.pressure(:, :, timeIndex).';
if options.Normalize
    scale = max(abs(field.pressure), [], "all");
    if scale > 0
        snapshot = snapshot / scale;
    end
end

if isempty(options.Parent)
    figure("Color", "w");
    ax = axes;
else
    ax = options.Parent;
    fig = ancestor(ax, "figure");
    if ~isempty(fig)
        set(fig, "Color", "w");
    end
end

imagesc(ax, field.r, field.z, snapshot);
set(ax, "YDir", "normal");
set(ax, "Color", "w", "XColor", "k", "YColor", "k", "FontSize", 11);
axis(ax, "image");
xlabel(ax, "radius r [m]");
ylabel(ax, "distance z [m]");
title(ax, sprintf("Step-struck drum field, t = %.3f ms", 1e3 * field.t(timeIndex)), ...
    "Color", "k");
cb = colorbar(ax);
set(cb, "Color", "k");
if options.Normalize
    ylabel(cb, "normalized pressure");
else
    ylabel(cb, "pressure [Pa]");
end
colormap(ax, "turbo");
end
