function mesh = readVolTriTet(volFile)
%READVOLTRITET Read the tri/tet subset of a Netgen .vol file.
%
% The MATLAB/Lukas/Gypsilab lane deliberately accepts only boundary triangles
% and volume tetrahedra. Quad/hex/prism/pyramid records are rejected instead of
% being split implicitly.

arguments
    volFile (1,1) string
end

txt = fileread(volFile);
lines = strip(splitlines(txt));

mesh = struct();
mesh.source = volFile;
mesh.vtx = zeros(0, 3);
mesh.tri = zeros(0, 3);
mesh.triCol = zeros(0, 1);
mesh.tet = zeros(0, 4);
mesh.tetMat = zeros(0, 1);
mesh.materials = containers.Map('KeyType', 'double', 'ValueType', 'char');
mesh.boundaries = containers.Map('KeyType', 'double', 'ValueType', 'char');
mesh.policy = "netgen_vol_tri_tet_only_shared_one_based_nodes";

i = 1;
while i <= numel(lines)
    key = nextDataLine();
    if key == ""
        break
    end

    switch char(lower(key))
        case 'mesh3d'
            % Header only.

        case 'endmesh'
            break

        case 'dimension'
            dim = str2double(nextDataLine());
            if dim ~= 3
                error("readVolTriTet:dimension", ".vol dimension must be 3.");
            end

        case 'geomtype'
            nextDataLine();

        case 'facedescriptors'
            skipCountedSection();

        case 'surfaceelements'
            n = readCount("surfaceelements");
            tri = zeros(n, 3);
            col = zeros(n, 1);
            for k = 1:n
                vals = sscanf(char(nextDataLine()), "%d").';
                if numel(vals) < 5
                    error("readVolTriTet:surface", "Surface record %d is too short.", k);
                end
                nodeCount = vals(5);
                if nodeCount ~= 3
                    error("readVolTriTet:surface", ...
                        "Only triangle surface elements are accepted; record %d has %d nodes.", k, nodeCount);
                end
                if numel(vals) < 5 + nodeCount
                    error("readVolTriTet:surface", "Surface record %d node count is inconsistent.", k);
                end
                col(k) = vals(2);
                tri(k, :) = vals(6:8);
            end
            mesh.tri = tri;
            mesh.triCol = col;

        case 'volumeelements'
            n = readCount("volumeelements");
            tet = zeros(n, 4);
            mat = zeros(n, 1);
            for k = 1:n
                vals = sscanf(char(nextDataLine()), "%d").';
                if numel(vals) < 2
                    error("readVolTriTet:volume", "Volume record %d is too short.", k);
                end
                nodeCount = vals(2);
                if nodeCount ~= 4
                    error("readVolTriTet:volume", ...
                        "Only tetrahedral volume elements are accepted; record %d has %d nodes.", k, nodeCount);
                end
                if numel(vals) < 2 + nodeCount
                    error("readVolTriTet:volume", "Volume record %d node count is inconsistent.", k);
                end
                mat(k) = vals(1);
                tet(k, :) = vals(3:6);
            end
            mesh.tet = tet;
            mesh.tetMat = mat;

        case 'points'
            n = readCount("points");
            vtx = zeros(n, 3);
            for k = 1:n
                vals = sscanf(char(nextDataLine()), "%f").';
                if numel(vals) < 3
                    error("readVolTriTet:points", "Point record %d must contain x y z.", k);
                end
                vtx(k, :) = vals(1:3);
            end
            mesh.vtx = vtx;

        case 'materials'
            n = readCount("materials");
            for k = 1:n
                parts = split(nextDataLine(), " ", 2);
                keyNum = str2double(parts(1));
                if isscalar(parts)
                    mesh.materials(keyNum) = sprintf("material_%d", keyNum);
                else
                    mesh.materials(keyNum) = char(parts(2));
                end
            end

        case 'bcnames'
            n = readCount("bcnames");
            for k = 1:n
                parts = split(nextDataLine(), " ", 2);
                keyNum = str2double(parts(1));
                if isscalar(parts)
                    mesh.boundaries(keyNum) = sprintf("boundary_%d", keyNum);
                else
                    mesh.boundaries(keyNum) = char(parts(2));
                end
            end

        case {'edgesegmentsgi2', 'pointelements', 'face_colours', 'face_transparencies'}
            skipCountedSection();

        case 'curvedelements'
            % High-order data may follow. The first-order tri/tet topology is
            % already sufficient for the first Lukas/Gypsilab coupling tests.
            break

        otherwise
            error("readVolTriTet:section", "Unsupported .vol section: %s", key);
    end
end

assertNodeReferences(mesh);
mesh.traceNodeIds = unique(mesh.tri(:));
mesh.summary = struct( ...
    "points", size(mesh.vtx, 1), ...
    "triangles", size(mesh.tri, 1), ...
    "tets", size(mesh.tet, 1), ...
    "materials", mesh.materials.Count, ...
    "boundaries", mesh.boundaries.Count);

    function line = nextDataLine()
        line = "";
        while i <= numel(lines)
            candidate = lines(i);
            i = i + 1;
            if candidate == "" || startsWith(candidate, "#")
                continue
            end
            line = candidate;
            return
        end
    end

    function n = readCount(context)
        raw = nextDataLine();
        n = str2double(raw);
        if isnan(n)
            error("readVolTriTet:count", "Invalid count for %s: %s", context, raw);
        end
    end

    function skipCountedSection()
        nSkip = readCount("counted section");
        for s = 1:nSkip
            nextDataLine();
        end
    end
end


function assertNodeReferences(mesh)
%ASSERTNODEREFERENCES Ensure all tri/tet connectivity points at .vol nodes.

n = size(mesh.vtx, 1);
allIds = [mesh.tri(:); mesh.tet(:)];
if any(allIds < 1) || any(allIds > n)
    error("readVolTriTet:nodes", "Connectivity references nodes outside 1..N.");
end
end
