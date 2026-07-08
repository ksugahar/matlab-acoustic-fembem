function tests = testCurvedPanelSphereConvergence
%TESTCURVEDPANELSPHERECONVERGENCE Curved vs flat panel BEM on the unit sphere.
%
% The straight-panel BEM's analytic deviation is faceted-geometry dominated
% (testHelmholtzScattering).  Replacing the flat panels by quadratic
% isoparametric curved panels -- edge nodes projected onto the sphere, the
% SOLUTION still P1 -- removes that O(h^2) faceting error: at the SAME mesh the
% curved lane is ~10-200x closer to the analytic truth.  This A/B locks the
% effect (Laplace capacitance -> 4*pi*R and sound-soft plane-wave scattering
% vs the partial-wave series) and writes a committed result JSON.
%
% Measured 2026-07-08 (R2026a), max curved errors:
%   capacitance  coarse 1.7e-4  fine 2.9e-5   (flat 1.3e-2 / 5.4e-3)
%   scatter k0.5 coarse 6.0e-4  fine 1.5e-4   (flat 1.4e-2 / 6.1e-3)
%   scatter k2.0 coarse 5.0e-3  fine 1.2e-3   (flat 3.6e-2 / 1.4e-2)
% curved/flat improvement 7x (k2 coarse) .. 190x (capacitance fine).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testCurvedBeatsFlatAndWritesJson(testCase)
R = 1.0;
proj = CurvedPanelQuadrature.sphereProjection(R);
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
meshes = ["unit_sphere_coarse", "unit_sphere_fine"];

% curved-error golden bands (measured x ~1.7 margin)
capBand = struct("unit_sphere_coarse", 5e-4, "unit_sphere_fine", 1e-4);
scatBand = struct( ...
    "unit_sphere_coarse", struct("k0p5", 1.2e-3, "k2", 8e-3), ...
    "unit_sphere_fine",   struct("k0p5", 4e-4,   "k2", 2.5e-3));

cases = struct("mesh", {}, "quantity", {}, "wavenumber", {}, ...
    "n_tri", {}, "n_node", {}, "flat_error", {}, "curved_error", {}, "improvement", {});

for name = meshes
    surface = fixtureSurface(name + ".vol");
    nT = size(surface.tri, 1); nN = size(surface.vtx, 1);

    % --- Laplace capacitance: potential 1 -> total charge 4*pi*R ---
    g0 = ones(nN, 1);
    solF = curvedSingleLayerDirichletSolve(surface, g0);
    solC = curvedSingleLayerDirichletSolve(surface, g0, "Projection", proj);
    verifyEqual(testCase, solF.status, "ok");
    verifyEqual(testCase, solC.status, "ok");
    eFcap = abs(solF.totalCharge - 4*pi*R) / (4*pi*R);
    eCcap = abs(solC.totalCharge - 4*pi*R) / (4*pi*R);
    verifyLessThan(testCase, eCcap, capBand.(name));
    verifyLessThan(testCase, eCcap, eFcap);
    verifyGreaterThan(testCase, eFcap / eCcap, 5);
    cases(end+1) = makeCase(name, "laplace_capacitance", 0.0, nT, nN, eFcap, eCcap); %#ok<AGROW>

    % --- sound-soft plane-wave scattering vs analytic series ---
    for k = [0.5 2.0]
        g = -exp(1i * k * surface.vtx(:, 3));
        solFk = curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k);
        solCk = curvedSingleLayerDirichletSolve(surface, g, "Wavenumber", k, ...
            "Projection", proj);
        ref = softSphereScattering(k, R, probes);
        verifyLessThan(testCase, ref.truncationTail, 1e-12);
        eFk = max(abs(solFk.potentialAt(probes) - ref.scattered) ./ abs(ref.scattered));
        eCk = max(abs(solCk.potentialAt(probes) - ref.scattered) ./ abs(ref.scattered));
        band = scatBand.(name).("k" + strrep(string(k), ".", "p"));
        verifyLessThan(testCase, eCk, band);
        verifyLessThan(testCase, eCk, eFk);
        verifyGreaterThan(testCase, eFk / eCk, 5);
        cases(end+1) = makeCase(name, "soft_sphere_scatter", k, nT, nN, eFk, eCk); %#ok<AGROW>
    end
end

% --- refinement makes the curved lane monotonically better ---
for q = ["laplace_capacitance", "soft_sphere_scatter"]
    for k = unique([cases.wavenumber])
        ec = arrayfun(@(c) c.curved_error, ...
            cases(strcmp({cases.quantity}, q) & [cases.wavenumber] == k));
        if numel(ec) == 2
            verifyLessThan(testCase, ec(2), ec(1));   % fine < coarse
        end
    end
end

writeResultJson(cases, R, probes);
end


function c = makeCase(name, quantity, k, nT, nN, eF, eC)
c = struct("mesh", string(name), "quantity", string(quantity), ...
    "wavenumber", k, "n_tri", nT, "n_node", nN, ...
    "flat_error", eF, "curved_error", eC, "improvement", eF / eC);
end


function writeResultJson(cases, R, probes)
result = struct();
result.kind = "curved_panel_sphere_convergence";
result.policy = "isoparametric_curved_panel_removes_flat_faceting_error";
result.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
result.matlab_version = string(version);
result.hostname = string(getenv("COMPUTERNAME"));
result.radius = R;
result.probe_points = probes;
result.quadrature_order = 7;
result.duffy_order = 6;
result.cases = cases;
here = fileparts(mfilename("fullpath"));
fid = fopen(fullfile(here, "curvedPanelSphereConvergence.json"), "w");
fwrite(fid, jsonencode(result, "PrettyPrint", true));
fclose(fid);
end


function surface = fixtureSurface(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
file = fullfile(repoRoot, "fixtures", "mesh_topology", name);
surface = VolMesh(file).boundary();
end
