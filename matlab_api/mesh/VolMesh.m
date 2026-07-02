classdef VolMesh
%VOLMESH First-order tri/tet volume mesh read from a Netgen .vol file.
%
% The mesh is the parsed .vol content plus its source identity:
%
%   mesh = VolMesh("model.vol");
%   mesh.vtx           % node coordinates (nNodes x 3)
%   mesh.tet           % tetrahedra as one-based volume node ids (nTets x 4)
%   mesh.tri           % boundary triangles as volume node ids (nTris x 3)
%   mesh.boundary()    % compact SurfaceMesh for the BEM spaces
%
% Only triangle surfaces and tetrahedral volumes are accepted; the parser
% readVolTriTet rejects quad/hex/prism/pyramid and curved records instead of
% splitting them silently.

properties
    vtx                 % node coordinates (nNodes x 3)
    tet                 % tetrahedra (nTets x 4), one-based volume node ids
    tetMat              % Netgen material number per tetrahedron (nTets x 1)
    tri                 % boundary triangles (nTris x 3), one-based volume node ids
    triCol              % Netgen boundary number per triangle (nTris x 1)
    materials           % containers.Map: material number -> name
    boundaries          % containers.Map: boundary number -> name
    boundaryOrientation % stored-normal vs outward classification report
    traceNodeIds        % sorted volume node ids that lie on the boundary
    policy              % mesh intake policy id
    summary             % point/triangle/tet/label counts
    sourcePath          % .vol file path
    sourceFileId        % "sha256:..." identity of the .vol bytes
    meshId              % "netgen_vol:<file>" mesh identity
end

methods
    function mesh = VolMesh(volFile)
        arguments
            volFile (1,1) string
        end
        raw = readVolTriTet(volFile);
        mesh.vtx = raw.vtx;
        mesh.tet = raw.tet;
        mesh.tetMat = raw.tetMat;
        mesh.tri = raw.tri;
        mesh.triCol = raw.triCol;
        mesh.materials = raw.materials;
        mesh.boundaries = raw.boundaries;
        mesh.boundaryOrientation = raw.boundaryOrientation;
        mesh.traceNodeIds = raw.traceNodeIds;
        mesh.policy = raw.policy;
        mesh.summary = raw.summary;
        mesh.sourcePath = string(volFile);
        mesh.sourceFileId = volFileSha256Id(volFile);
        [~, stem, ext] = fileparts(mesh.sourcePath);
        mesh.meshId = "netgen_vol:" + string(stem) + string(ext);
    end

    function surface = boundary(mesh)
        %BOUNDARY Compact boundary triangle mesh for the BEM spaces.
        surface = SurfaceMesh(mesh);
    end

    function names = boundaryNames(mesh)
        %BOUNDARYNAMES Expand Netgen boundary numbers to per-triangle labels.
        names = strings(numel(mesh.triCol), 1);
        for k = 1:numel(mesh.triCol)
            id = mesh.triCol(k);
            if isKey(mesh.boundaries, id)
                names(k) = string(mesh.boundaries(id));
            else
                names(k) = "boundary_" + id;
            end
        end
    end

    function r = regions(mesh)
        %REGIONS Group tetrahedra by Netgen material number.
        matIds = unique(mesh.tetMat(:)).';
        r = struct("matnr", {}, "name", {}, "Elements", {}, "Nodes", {});
        for k = 1:numel(matIds)
            matnr = matIds(k);
            elems = find(mesh.tetMat == matnr);
            nodes = unique(mesh.tet(elems, :));
            r(k).matnr = matnr;
            if isKey(mesh.materials, matnr)
                r(k).name = string(mesh.materials(matnr));
            else
                r(k).name = "material_" + matnr;
            end
            r(k).Elements = elems(:);
            r(k).Nodes = nodes(:);
        end
    end
end
end


function id = volFileSha256Id(path)
%VOLFILESHA256ID Stable source identity for the .vol-to-trace handoff.

fid = fopen(path, "rb");
if fid < 0
    error("VolMesh:file", "Cannot open .vol source: %s", path);
end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, "*uint8");
md = javaMethod("getInstance", "java.security.MessageDigest", "SHA-256");
md.update(typecast(bytes(:), "int8"));
hash = typecast(md.digest(), "uint8");
hex = lower(reshape(dec2hex(hash, 2).', 1, []));
id = "sha256:" + string(hex);
clear cleanup
end
