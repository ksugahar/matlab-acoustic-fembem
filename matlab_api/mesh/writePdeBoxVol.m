function report = writePdeBoxVol(volFile, options)
%WRITEPDEBOXVOL Create a simple PDE Toolbox box mesh and export it to .vol.
%
%   report = writePdeBoxVol("box.vol", Hmax=0.25)
%
% A box convenience wrapper over writePdeGeometryVol (multicuboid geometry).
% It needs the optional MATLAB PDE Toolbox; writePdeMeshVol (an existing mesh)
% and structuredBoxVol (a box, no toolbox) are the toolbox-free alternatives.

arguments
    volFile (1,1) string
    options.Size (1,3) double {mustBePositive} = [1 1 1]
    options.Hmax (1,1) double {mustBePositive} = 0.25
    options.MaterialName (1,1) string = "domain"
    options.BoundaryName (1,1) string = "outer"
end

% Check availability without the "file" filter: in R2026a multicuboid is a
% BUILT-IN (exist == 5), so the old exist(...,"file")==2 test wrongly rejected an
% installed PDE Toolbox.  exist(name)==0 means the function is absent entirely.
if exist("createpde") == 0 || exist("multicuboid") == 0
    error("writePdeBoxVol:pdeToolboxUnavailable", ...
        "PDE Toolbox is required for writePdeBoxVol. Use writePdeMeshVol with an existing mesh instead.");
end

report = writePdeGeometryVol(volFile, ...
    multicuboid(options.Size(1), options.Size(2), options.Size(3)), ...
    Hmax=options.Hmax, MaterialName=options.MaterialName, BoundaryName=options.BoundaryName);
report.generator = "pde_toolbox_multicuboid";
report.box_size = options.Size;
end
