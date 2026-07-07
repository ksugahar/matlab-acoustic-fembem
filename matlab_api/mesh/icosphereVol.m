function report = icosphereVol(volFile, options)
%ICOSPHEREVOL Write a sphere .vol by subdividing an icosahedron (pure MATLAB).
%
%   icosphereVol("sphere.vol")                       % unit sphere, 3 subdivisions
%   icosphereVol("s.vol", Subdivisions=4)            % finer (5120 tris, 2562 nodes)
%   icosphereVol("s.vol", Radius=2, Center=[0 0 1])  % moved / scaled
%   icosphereVol("s.vol", SurfaceOnly=true)          % boundary triangles only
%
% A pure-MATLAB sphere mesher: subdivide the icosahedron, project each vertex
% to the sphere, orient every triangle outward, and (by default) add one
% interior node at the center with a "cone" tetrahedron per face so the .vol
% is a closed, outward-oriented volume mesh.  No PDE Toolbox, no netgen --
% MATLAB alone generates the mesh and writeVol writes it.  Subdivisions s give
% 20*4^s triangles and 10*4^s+2 surface nodes (s=3 -> 1280 tris, 642 nodes).

arguments
    volFile (1,1) string
    options.Radius (1,1) double {mustBePositive} = 1.0
    options.Subdivisions (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    options.Center (1,3) double = [0 0 0]
    options.BoundaryName (1,1) string = "outer"
    options.MaterialName (1,1) string = "domain"
    options.SurfaceOnly (1,1) logical = false
end

[V, F] = icosahedron();
for s = 1:options.Subdivisions
    [V, F] = subdivideTriangles(V, F);
end
V = V ./ vecnorm(V, 2, 2);                          % project vertices to the unit sphere

% orient every face outward (normal points away from the origin)
faceNormal = cross(V(F(:, 2), :) - V(F(:, 1), :), V(F(:, 3), :) - V(F(:, 1), :), 2);
faceCentroid = (V(F(:, 1), :) + V(F(:, 2), :) + V(F(:, 3), :)) / 3;
inward = sum(faceNormal .* faceCentroid, 2) < 0;
F(inward, [2 3]) = F(inward, [3 2]);

pts = options.Center + options.Radius * V;

if options.SurfaceOnly
    report = writeVol(volFile, pts, F, ...
        BoundaryNames=options.BoundaryName, MaterialNames=options.MaterialName);
else
    centerId = size(pts, 1) + 1;
    pts = [pts; options.Center];                    % one interior node
    tets = [F, repmat(centerId, size(F, 1), 1)];    % a cone tet under each face
    report = writeVol(volFile, pts, F, Tets=tets, ...
        BoundaryNames=options.BoundaryName, MaterialNames=options.MaterialName);
end

report.tool = "icosphere_vol";
report.subdivisions = options.Subdivisions;
report.radius = options.Radius;
report.surface_nodes = size(V, 1);
end


function [V, F] = icosahedron()
%ICOSAHEDRON 12 vertices, 20 faces of a regular icosahedron.
t = (1 + sqrt(5)) / 2;
V = [-1  t  0;  1  t  0; -1 -t  0;  1 -t  0;
      0 -1  t;  0  1  t;  0 -1 -t;  0  1 -t;
      t  0 -1;  t  0  1; -t  0 -1; -t  0  1];
F = [ 1 12  6;  1  6  2;  1  2  8;  1  8 11;  1 11 12;
      2  6 10;  6 12  5; 12 11  3; 11  8  7;  8  2  9;
      4 10  5;  4  5  3;  4  3  7;  4  7  9;  4  9 10;
      5 10  6;  3  5 12;  7  3 11;  9  7  8; 10  9  2];
end


function [V, F] = subdivideTriangles(V, F)
%SUBDIVIDETRIANGLES 1-to-4 triangle split with shared edge midpoints.
edgeMap = containers.Map("KeyType", "char", "ValueType", "double");
Vcell = V;
Fnew = zeros(4 * size(F, 1), 3);
row = 0;
for k = 1:size(F, 1)
    a = F(k, 1); b = F(k, 2); c = F(k, 3);
    ab = midpoint(a, b); bc = midpoint(b, c); ca = midpoint(c, a);
    Fnew(row + 1, :) = [a ab ca];
    Fnew(row + 2, :) = [b bc ab];
    Fnew(row + 3, :) = [c ca bc];
    Fnew(row + 4, :) = [ab bc ca];
    row = row + 4;
end
V = Vcell;
F = Fnew;

    function m = midpoint(i, j)
        key = sprintf("%d_%d", min(i, j), max(i, j));
        if isKey(edgeMap, key)
            m = edgeMap(key);
            return
        end
        Vcell(end + 1, :) = 0.5 * (Vcell(i, :) + Vcell(j, :)); %#ok<AGROW>
        m = size(Vcell, 1);
        edgeMap(key) = m;
    end
end
