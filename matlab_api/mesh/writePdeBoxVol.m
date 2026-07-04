function report = writePdeBoxVol(volFile, options)
%WRITEPDEBOXVOL Create a simple PDE Toolbox box mesh and export it to .vol.
%
%   report = writePdeBoxVol("box.vol", Hmax=0.25)
%
% This convenience function is optional at runtime: it needs MATLAB PDE
% Toolbox.  The core exporter writePdeMeshVol can be tested without PDE
% Toolbox by passing a struct with Nodes and Elements fields.

arguments
    volFile (1,1) string
    options.Size (1,3) double {mustBePositive} = [1 1 1]
    options.Hmax (1,1) double {mustBePositive} = 0.25
    options.MaterialName (1,1) string = "domain"
    options.BoundaryName (1,1) string = "outer"
end

if exist("createpde", "file") ~= 2 || exist("multicuboid", "file") ~= 2
    error("writePdeBoxVol:pdeToolboxUnavailable", ...
        "PDE Toolbox is required for writePdeBoxVol. Use writePdeMeshVol with an existing mesh instead.");
end

model = createpde();
model.Geometry = multicuboid(options.Size(1), options.Size(2), options.Size(3));
generateMesh(model, "Hmax", options.Hmax, "GeometricOrder", "linear");

report = writePdeMeshVol(model.Mesh, volFile, ...
    MaterialName=options.MaterialName, ...
    BoundaryName=options.BoundaryName);
report.generator = "pde_toolbox_multicuboid";
report.box_size = options.Size;
report.hmax = options.Hmax;
end
