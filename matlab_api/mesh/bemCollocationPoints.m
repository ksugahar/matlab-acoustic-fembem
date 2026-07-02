function [points, weights] = bemCollocationPoints(input)
%BEMCOLLOCATIONPOINTS Collocation points and area weights for BEM operators.
%
% Accepts the three shapes the teaching operators work on:
%
%   Nx3 double     -> the points themselves, unit weights
%   SurfaceMesh    -> triangle centroids, triangle areas
%   FemBemModel    -> the centroids/areas of its boundary surface
%
% Anything else is rejected loudly; there is no silent conversion.

if isnumeric(input)
    points = input;
    weights = ones(size(points, 1), 1);
elseif isa(input, "SurfaceMesh")
    points = input.centroids();
    weights = input.areas();
elseif isa(input, "FemBemModel")
    points = input.surface.centroids();
    weights = input.surface.areas();
else
    error("bemCollocationPoints:input", ...
        "Input must be an Nx3 point array, a SurfaceMesh, or a FemBemModel.");
end
end
