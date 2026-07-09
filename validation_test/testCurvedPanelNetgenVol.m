function tests = testCurvedPanelNetgenVol
%TESTCURVEDPANELNETGENVOL Optional netgen path: curved .vol -> MATLAB high accuracy.
%
% NGSolve (which owns the curved-element engine) evaluates each boundary
% triangle's curved geometry at the Lagrange nodes and writes the physical
% COORDINATES to a convention-free JSON companion (the saved .vol stays P1;
% MATLAB never touches Netgen's curvedelements coefficient basis).  MATLAB reads
% both with curvedQuadratureFromNetgen and gets the SAME faceting-free accuracy
% as the self-generated analytic-projection lane -- but for GENERAL geometry.
%
% Fixture: fixtures/curved_panels/ng_sphere.vol (+ ng_sphere_curved_p2.json),
% a netgen OCC sphere curved to order 2 by tools/export_curved_boundary_nodes.py.
% Measured 2026-07-08 (R2026a): netgen-fed geomNodes match analytic projection to
% 5.5e-4; capacitance 3.4e-4, scatter k0.5 8.8e-4, k2 6.5e-3 -- ~40x/18x/6x below
% the flat lane and within ~1.6x of the analytic-projection lane.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testNetgenFedMatchesAnalyticAndBeatsFlatWritesJson(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
vol = fullfile(repoRoot, "fixtures", "curved_panels", "ng_sphere.vol");
json = fullfile(repoRoot, "fixtures", "curved_panels", "ng_sphere_curved_p2.json");
surface = VolMesh(vol).boundary();
R = mean(vecnorm(surface.vtx, 2, 2));
probes = [2 0 0; 0 0 3; -1.2 1.6 0];

proj = CurvedPanelQuadrature.sphereProjection(R);
qN = curvedQuadratureFromNetgen(surface, json, 7);
qA = CurvedPanelQuadrature(surface, 7, proj, 2);

% netgen curved nodes match the analytic projection to the mesh-curving accuracy
nodeDiff = max(abs(qN.geomNodes(:) - qA.geomNodes(:)));
verifyLessThan(testCase, nodeDiff, 2e-3);

% area: netgen-fed close to 4*pi*R^2 and far better than flat
areaTrue = 4 * pi * R^2;
areaFlat = sum(CurvedPanelQuadrature(surface, 7, @(X) X, 1).weights);
areaN = sum(qN.weights);
verifyLessThan(testCase, abs(areaN - areaTrue) / areaTrue, 5e-3);
verifyLessThan(testCase, abs(areaN - areaTrue), 0.1 * abs(areaFlat - areaTrue));

% capacitance + scattering: netgen-fed beats flat and tracks analytic projection
g0 = ones(size(surface.vtx, 1), 1);
capT = 4 * pi * R;
eFcap = relerr(curvedSingleLayerDirichletSolve(surface, g0, "CurveOrder", 1, "Projection", proj).totalCharge, capT);
eNcap = relerr(curvedSingleLayerDirichletSolve(surface, g0, "GeomNodes", qN.geomNodes).totalCharge, capT);
eAcap = relerr(curvedSingleLayerDirichletSolve(surface, g0, "Projection", proj, "CurveOrder", 2).totalCharge, capT);
verifyLessThan(testCase, eNcap, 1e-3);
verifyLessThan(testCase, eNcap, eFcap / 10);
verifyLessThan(testCase, eNcap, 3 * eAcap);

cases = struct("quantity", {"laplace_capacitance"}, "wavenumber", {0.0}, ...
    "flat_error", {eFcap}, "netgen_error", {eNcap}, "analytic_error", {eAcap});

scatBand = struct("k0p5", 2e-3, "k2", 1e-2);
for k = [0.5 2.0]
    g = -exp(1i * k * surface.vtx(:, 3));
    ref = softSphereScattering(k, R, probes).scattered;
    eF = relerrVec(curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k, "CurveOrder", 1, "Projection", proj).potentialAt(probes), ref);
    eN = relerrVec(curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k, "GeomNodes", qN.geomNodes).potentialAt(probes), ref);
    eA = relerrVec(curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k, "Projection", proj, "CurveOrder", 2).potentialAt(probes), ref);
    verifyLessThan(testCase, eN, scatBand.("k" + strrep(string(k), ".", "p")));
    verifyLessThan(testCase, eN, eF / 5);
    verifyLessThan(testCase, eN, 3 * eA);
    cases(end+1) = struct("quantity", "soft_sphere_scatter", "wavenumber", k, ...
        "flat_error", eF, "netgen_error", eN, "analytic_error", eA); %#ok<AGROW>
end

writeResultJson(cases, R, nodeDiff, size(surface.tri, 1), size(surface.vtx, 1));
end


function e = relerr(v, ref)
e = abs(v - ref) / abs(ref);
end


function e = relerrVec(v, ref)
e = max(abs(v - ref) ./ abs(ref));
end


function writeResultJson(cases, R, nodeDiff, nTri, nNode)
result = struct();
result.kind = "curved_panel_netgen_vol";
result.policy = "netgen_curved_mesh_to_matlab_via_convention_free_node_companion";
result.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
result.matlab_version = string(version);
result.hostname = string(getenv("COMPUTERNAME"));
result.fixture = "fixtures/curved_panels/ng_sphere.vol + ng_sphere_curved_p2.json";
result.radius = R;
result.n_tri = nTri;
result.n_node = nNode;
result.netgen_vs_analytic_node_maxdiff = nodeDiff;
result.cases = cases;
here = fileparts(mfilename("fullpath"));
fid = fopen(fullfile(here, "curvedPanelNetgenVol.json"), "w");
fwrite(fid, jsonencode(result, "PrettyPrint", true));
fclose(fid);
end
