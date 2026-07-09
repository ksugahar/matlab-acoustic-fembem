function tests = testCurvedPanelCurveOrderSweep
%TESTCURVEDPANELCURVEORDERSWEEP Curve order is the geometry lever (fes fixed at P1).
%
% Sweeping the CURVE order 1 (flat) -> 2 (quadratic) -> 3 (cubic) with the
% SOLUTION fixed at P1 shows the two error channels explicitly:
%   - Laplace capacitance (smooth): keeps improving with curve order (1.3e-2 ->
%     1.7e-4 -> 4.7e-5 coarse) -- geometry is the bottleneck the whole way, so
%     curve order is the lever.
%   - Sound-soft scattering k=2: 1->2 is a big jump (7x) but 2->3 PLATEAUS
%     (5.0e-3 -> 4.7e-3) -- once the geometry is resolved the P1 (fes) density
%     error is the floor, and curve order stops helping.
% Together: curve order dominates while geometry is the binding constraint;
% after that, fes order takes over.  Measured 2026-07-08 (R2026a).

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testCurveOrderSweepAndWritesJson(testCase)
proj = CurvedPanelQuadrature.sphereProjection(1.0);
probes = [2 0 0; 0 0 3; -1.2 1.6 0];
kScatter = 2.0;

cap = struct("unit_sphere_coarse", [], "unit_sphere_fine", []);
sca = struct("unit_sphere_coarse", [], "unit_sphere_fine", []);
cases = struct("mesh", {}, "quantity", {}, "curve_order", {}, ...
    "n_tri", {}, "error", {});

for name = ["unit_sphere_coarse", "unit_sphere_fine"]
    surface = fixtureSurface(name + ".vol");
    nT = size(surface.tri, 1);
    g0 = ones(size(surface.vtx, 1), 1);
    gS = -exp(1i * kScatter * surface.vtx(:, 3));
    ref = softSphereScattering(kScatter, 1.0, probes);
    ecap = zeros(1, 3); esca = zeros(1, 3);
    for co = [1 2 3]
        solC = curvedSingleLayerDirichletSolve(surface, g0, "Projection", proj, "CurveOrder", co);
        ecap(co) = abs(solC.totalCharge - 4*pi) / (4*pi);
        solS = curvedSingleLayerDirichletSolve(surface, gS, "Wavenumber", kScatter, ...
            "Projection", proj, "CurveOrder", co);
        esca(co) = max(abs(solS.potentialAt(probes) - ref.scattered) ./ abs(ref.scattered));
        cases(end+1) = struct("mesh", string(name), "quantity", "laplace_capacitance", ...
            "curve_order", co, "n_tri", nT, "error", ecap(co)); %#ok<AGROW>
        cases(end+1) = struct("mesh", string(name), "quantity", "soft_sphere_scatter", ...
            "curve_order", co, "n_tri", nT, "error", esca(co)); %#ok<AGROW>
    end
    cap.(name) = ecap; sca.(name) = esca;

    % capacitance: geometry-limited -> every curve-order bump keeps improving
    verifyLessThan(testCase, ecap(2), ecap(1) / 5);
    verifyLessThan(testCase, ecap(3), ecap(2) / 2);

    % scattering: 1->2 big, then 2->3 PLATEAUS (P1 fes floor reached)
    verifyLessThan(testCase, esca(2), esca(1) / 5);   % curving pays off
    verifyLessThan(testCase, esca(3), esca(2) * 1.05); % never worse
    verifyGreaterThan(testCase, esca(3), esca(2) / 1.5); % ...but plateaus
end

writeResultJson(cases, kScatter);
end


function writeResultJson(cases, kScatter)
result = struct();
result.kind = "curved_panel_curve_order_sweep";
result.policy = "curve_order_is_the_geometry_lever_until_fes_floor";
result.lesson = "capacitance keeps improving with curve order (geometry-limited); " + ...
    "k=2 scattering plateaus at curve order 2 (P1 fes density floor reached).";
result.generated_at_utc = string(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'"));
result.matlab_version = string(version);
result.hostname = string(getenv("COMPUTERNAME"));
result.radius = 1.0;
result.scatter_wavenumber = kScatter;
result.fes_order = 1;
result.cases = cases;
here = fileparts(mfilename("fullpath"));
fid = fopen(fullfile(here, "curvedPanelCurveOrderSweep.json"), "w");
fwrite(fid, jsonencode(result, "PrettyPrint", true));
fclose(fid);
end


function surface = fixtureSurface(name)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
file = fullfile(repoRoot, "fixtures", "mesh_topology", name);
surface = VolMesh(file).boundary();
end
