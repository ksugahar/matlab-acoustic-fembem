function [points, weights] = bemCollocationPoints(input)
%BEMCOLLOCATIONPOINTS P1 nodal points and lumped area weights for BEM.
%
% Accepts the three shapes the teaching operators work on:
%
%   Nx3 double     -> the points themselves, unit weights
%   SurfaceMesh    -> compact boundary nodes, lumped P1 nodal areas
%   FemBemModel    -> the P1 nodes/weights of its boundary surface
%
% Anything else is rejected loudly; there is no silent conversion.

if isnumeric(input)
    points = input;
    weights = ones(size(points, 1), 1);
elseif isa(input, "SurfaceMesh")
    points = input.vtx;
    weights = surfaceP1LumpedWeights(input);
elseif isa(input, "FemBemModel")
    points = input.surface.vtx;
    weights = surfaceP1LumpedWeights(input.surface);
else
    error("bemCollocationPoints:input", ...
        "Input must be an Nx3 point array, a SurfaceMesh, or a FemBemModel.");
end
end


function weights = surfaceP1LumpedWeights(surface)
%SURFACEP1LUMPEDWEIGHTS Integral of each P1 nodal basis over the surface.

nNodes = size(surface.vtx, 1);
weights = zeros(nNodes, 1);
areas = surface.areas();
for e = 1:size(surface.tri, 1)
    weights(surface.tri(e, :)) = weights(surface.tri(e, :)) + areas(e) / 3;
end
end
