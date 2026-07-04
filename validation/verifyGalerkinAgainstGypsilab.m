function report = verifyGalerkinAgainstGypsilab(volFile, gypsilabRoot)
%VERIFYGALERKINAGAINSTGYPSILAB Same-mesh check against the real Gypsilab.
%
%   report = verifyGalerkinAgainstGypsilab( ...
%       "fixtures/mesh_topology/unit_sphere_coarse.vol", ...
%       fullfile("path", "to", "gypsilab"));
%
% Builds the boundary P1 Galerkin single layer twice on the same mesh:
% GalerkinSingleLayer (this repo) versus Gypsilab's
% integral(sigma,sigma,u,G,v) + regularize(sigma,sigma,u,'[1/r]',v), then
% compares operator, mass, and the g = 1 capacitance. Requires the real
% Gypsilab source tree and errors loudly when it is missing (no fallback).

arguments
    volFile (1,1) string
    gypsilabRoot (1,1) string
end

for folder = ["openMsh", "openDom", "openFem", "openOpr", "miscellaneous"]
    p = fullfile(gypsilabRoot, folder);
    if ~isfolder(p)
        error("verifyGalerkinAgainstGypsilab:root", ...
            "Gypsilab folder not found: %s", p);
    end
    addpath(p);
end

mesh = VolMesh(volFile);
surface = mesh.boundary();
nNodes = size(surface.vtx, 1);

ours = GalerkinSingleLayer(surface, "QuadratureOrder", 3);
solOurs = singleLayerDirichletSolve(surface, ones(nNodes, 1));
[Mours, ~] = SurfaceP1Space(surface).mass();

mgy = surface.gypsilabMsh();
if size(mgy.vtx, 1) ~= nNodes
    error("verifyGalerkinAgainstGypsilab:vertices", ...
        "Gypsilab msh changed the vertex count.");
end
[tf, perm] = ismember(round(surface.vtx, 12), round(mgy.vtx, 12), "rows");
if ~all(tf)
    error("verifyGalerkinAgainstGypsilab:vertices", ...
        "Vertex matching between SurfaceMesh and Gypsilab msh failed.");
end

sigma = dom(mgy, 3);
u = fem(mgy, 'P1');
v = fem(mgy, 'P1');
Gxy = @(X, Y) femGreenKernel(X, Y, '[1/r]', 0);
Vgy = 1/(4*pi) .* integral(sigma, sigma, u, Gxy, v) ...
    + 1/(4*pi) .* regularize(sigma, sigma, u, '[1/r]', v);
Mgy = integral(sigma, u, v);
Vg = full(Vgy(perm, perm));
Mg = full(Mgy(perm, perm));

qgy = Vg \ (Mg * ones(nNodes, 1));
capacitanceGypsilab = sum(Mg * qgy);

report = struct();
report.kind = "galerkin_single_layer_gypsilab_cross_check";
report.volFile = string(volFile);
report.nNodes = nNodes;
report.operatorRelDiff = norm(ours.matrix - Vg, "fro") / norm(Vg, "fro");
report.massRelDiff = norm(full(Mours) - Mg, "fro") / norm(Mg, "fro");
report.capacitanceOurs = solOurs.totalCharge;
report.capacitanceGypsilab = capacitanceGypsilab;
report.capacitanceRelDiff = abs(solOurs.totalCharge - capacitanceGypsilab) ...
    / abs(capacitanceGypsilab);
report.checks = struct( ...
    "massIdentical", report.massRelDiff < 1e-12, ...
    "operatorClose", report.operatorRelDiff < 5e-3, ...
    "capacitanceClose", report.capacitanceRelDiff < 1e-3);
if all(structfun(@(x) logical(x), report.checks))
    report.status = "ok";
else
    report.status = "needs_attention";
end
end
