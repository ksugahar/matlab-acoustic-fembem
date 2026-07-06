function report = writePdeGeometryVol(volFile, geometry, options)
%WRITEPDEGEOMETRYVOL Mesh any PDE Toolbox geometry and export it to .vol.
%
%   report = writePdeGeometryVol("sphere.vol", multisphere(0.5), Hmax=0.12)
%   report = writePdeGeometryVol("cyl.vol",    multicylinder(0.4, 1.0), Hmax=0.15)
%   report = writePdeGeometryVol("box.vol",    multicuboid(1.2, 1.0, 0.8), Hmax=0.35)
%   g = importGeometry("part.step");
%   report = writePdeGeometryVol("part.vol", g, Hmax=0.05)
%
% Creates a first-order tetrahedral mesh of ANY MATLAB PDE Toolbox geometry
% (multicuboid / multisphere / multicylinder / importGeometry result) and exports
% it to the first-order Netgen .vol contract via writePdeMeshVol -- the general
% PDE-Toolbox-to-.vol path for arbitrary scatterer shapes.  Needs the OPTIONAL
% PDE Toolbox at runtime; writePdeMeshVol itself (an existing mesh) does not, and
% structuredBoxVol produces a box .vol with no toolbox at all.

arguments
    volFile (1,1) string
    geometry
    options.Hmax (1,1) double {mustBePositive} = 0.25
    options.MaterialName (1,1) string = "domain"
    options.BoundaryName (1,1) string = "outer"
end

% Availability check without the "file" filter (PDE geometry builders are
% built-ins in R2026a, so exist(...,"file") would wrongly reject them).
if exist("createpde") == 0
    error("writePdeGeometryVol:pdeToolboxUnavailable", ...
        "PDE Toolbox is required for writePdeGeometryVol. Use writePdeMeshVol " + ...
        "with an existing mesh, or structuredBoxVol for a box, instead.");
end

model = createpde();
model.Geometry = geometry;
generateMesh(model, "Hmax", options.Hmax, "GeometricOrder", "linear");

report = writePdeMeshVol(model.Mesh, volFile, ...
    MaterialName=options.MaterialName, BoundaryName=options.BoundaryName);
report.generator = "pde_toolbox_geometry";
report.hmax = options.Hmax;
end
