function report = verifyGalerkinAgainstNgsolve(volFile, matFile)
%VERIFYGALERKINAGAINSTNGSOLVE Same-mesh check against NGSolve's ngsolve.bem.
%
%   report = verifyGalerkinAgainstNgsolve( ...
%       "fixtures/mesh_topology/unit_sphere_coarse.vol");
%
% Compares this repo's Galerkin boundary-P1 Laplace operators against dense
% reference matrices assembled by ngsolve.bem on the SAME .vol mesh and the
% SAME continuous-P1 space (exported once by exportNgsolveBemReference.py;
% NGSolve H1 order-1 dof numbering equals the .vol point order, so the only
% reindexing is SurfaceMesh.volNodeIds - verified on the vertex coordinates
% before anything else). Checked:
%
%   mass          identical surface P1 mass            (< 1e-12, meas 1e-16)
%   V  (gss = 7)  single layer vs Sauter-Schwab ref    (< 1e-3, meas 3-4e-4)
%   K  (gss = 7)  outward-normal PV double layer       (< 8e-3, meas 3.4e-3)
%   capacitance   g = 1 solve, our V vs ngbem V + same mass
%                                                      (< 1e-3, meas 2e-4)
%
% ngsolve.bem shares our operator conventions exactly (measured: K[1] = -1/2
% to 1e-9, K[Y_1] -> -1/6, principal value on the diagonal); the remaining
% difference is our one-sided test Gauss quadrature (V reldiff converges
% 7e-2 / 7e-3 / 3-4e-4 at gss 1/3/7 - so the 1.1e-4 agreement with the real
% Gypsilab at gss 3 was same-quadrature error cancellation, while this check
% measures the true assembly error against a converged reference).
% Errors loudly when the reference .mat is missing (no fallback).

arguments
    volFile (1,1) string
    matFile (1,1) string = ""
end

if matFile == ""
    [~, base] = fileparts(volFile);
    matFile = fullfile(fileparts(mfilename("fullpath")), "data", ...
        "ngbem_reference_" + base + ".mat");
end
if ~isfile(matFile)
    error("verifyGalerkinAgainstNgsolve:reference", ...
        "ngsolve.bem reference not found: %s\nRegenerate with:\n  python validation/exportNgsolveBemReference.py %s %s", ...
        matFile, volFile, matFile);
end
S = load(matFile);

mesh = VolMesh(volFile);
surface = mesh.boundary();
ids = surface.volNodeIds;
nNodes = size(surface.vtx, 1);

vertexMismatch = max(abs(surface.vtx - S.vtx(ids, :)), [], "all");
if vertexMismatch > 1e-12
    error("verifyGalerkinAgainstNgsolve:vertices", ...
        "NGSolve vertex order does not match the .vol point order (max %.3e).", ...
        vertexMismatch);
end

Vng = S.V(ids, ids);
Kng = S.K(ids, ids);
Mng = S.M(ids, ids);

[Mours, ~] = SurfaceP1Space(surface).mass();
Vours = GalerkinSingleLayer(surface, "QuadratureOrder", 7).matrix;
Kours = GalerkinDoubleLayer(surface, "QuadratureOrder", 7).matrix;

solOurs = singleLayerDirichletSolve(surface, ones(nNodes, 1));
qng = Vng \ (full(Mours) * ones(nNodes, 1));
capacitanceNgsolve = sum(full(Mours) * qng);

report = struct();
report.kind = "galerkin_laplace_ngsolve_bem_cross_check";
report.volFile = string(volFile);
report.matFile = string(matFile);
report.nNodes = nNodes;
report.ngsolveVersion = string(S.ngsolveVersion);
report.referenceIntorder = S.intorder;
report.referenceIntorderConvergenceV = S.intorderConvergenceV;
report.referenceIntorderConvergenceK = S.intorderConvergenceK;
report.quadratureOrderOperators = 7;
report.quadratureOrderSolve = solOurs.quadratureOrder;
report.massRelDiff = norm(full(Mours) - Mng, "fro") / norm(Mng, "fro");
report.operatorVRelDiff = norm(Vours - Vng, "fro") / norm(Vng, "fro");
report.operatorKRelDiff = norm(Kours - Kng, "fro") / norm(Kng, "fro");
report.capacitanceOurs = solOurs.totalCharge;
report.capacitanceNgsolve = capacitanceNgsolve;
report.capacitanceRelDiff = abs(solOurs.totalCharge - capacitanceNgsolve) ...
    / abs(capacitanceNgsolve);
report.checks = struct( ...
    "referenceConverged", ...
        S.intorderConvergenceV < 1e-6 && S.intorderConvergenceK < 1e-6, ...
    "massIdentical", report.massRelDiff < 1e-12, ...
    "singleLayerClose", report.operatorVRelDiff < 1e-3, ...
    "doubleLayerClose", report.operatorKRelDiff < 8e-3, ...
    "capacitanceClose", report.capacitanceRelDiff < 1e-3);
if all(structfun(@(x) logical(x), report.checks))
    report.status = "ok";
else
    report.status = "needs_attention";
end
end
