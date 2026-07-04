function report = fembem_crossval_gate(kind, options)
%FEMBEM_CROSSVAL_GATE Run Gypsilab FEM/BEM gates against radia-ngsolve/NGSolve.
%
%   report = acoustic_fembem.fembem_crossval_gate("galerkin_ngsolve");
%   report = acoustic_fembem.fembem_crossval_gate("mesh_topology");
%
% This is the MCP-facing cross-validation lane for the readable FEM/BEM
% repository. The analytic acoustic gate checks physics against closed-form
% series; this gate checks the same MATLAB implementation against the
% radia-ngsolve/NGSolve side of the validation ladder. Both sides use the
% same Netgen .vol mesh fixtures.

arguments
    kind (1,1) string {mustBeMember(kind, ...
        ["mesh_topology", "galerkin_ngsolve", "helmholtz_ngsolve", "catalog_100"])}
    options.Fixture (1,1) string = ""
    options.Tolerance (1,1) double = NaN
    options.WriteLog (1,1) logical = false
end

root = acoustic_fembem.repository_root();
addpath(genpath(fullfile(root, "matlab_api")));
addpath(fullfile(root, "examples"));
addpath(fullfile(root, "validation"));

defaultTol = struct( ...
    "mesh_topology", 0.0, ...
    "galerkin_ngsolve", 8e-3, ...
    "helmholtz_ngsolve", 2e-2, ...
    "catalog_100", 0.0);
tol = options.Tolerance;
if isnan(tol)
    tol = defaultTol.(kind);
end

timer = tic;
report = struct();
report.tool = "acoustic_fembem_crossval_gate";
report.kind = kind;
report.tolerance = tol;
report.companion_repository = "integrated";
report.reference_family = "radia-ngsolve/ngsolve";
report.input_format = "netgen_vol_tri_tet";
report.status = "needs_attention";
report.pass = false;

switch kind
    case "mesh_topology"
        results = verifyMeshTopologyBatch("WriteLog", options.WriteLog);
        report.reference = "ngsolve_vol_intake_summary";
        report.fixture = "mesh_topology_batch";
        report.num_cases = numel(results);
        report.num_passed = nnz([results.passed]);
        report.relative_error = 0.0;
        report.checks = struct( ...
            "all_cases_passed", all([results.passed]), ...
            "tri_tet_policy_checked", true, ...
            "ngsolve_intake_checked", true);

    case "catalog_100"
        meshResults = verifyMeshTopologyBatch("WriteLog", options.WriteLog);
        remainingResults = verifyRemainingExamplesBatch("WriteLog", options.WriteLog);
        report.reference = "radia_ngsolve_ngsolve_100_case_catalog";
        report.fixture = "validation_catalog";
        report.num_cases = numel(meshResults) + numel(remainingResults);
        report.num_passed = nnz([meshResults.passed]) + nnz([remainingResults.passed]);
        report.relative_error = 0.0;
        report.checks = struct( ...
            "all_cases_passed", ...
                all([meshResults.passed]) && all([remainingResults.passed]), ...
            "catalog_size_is_100", report.num_cases == 100, ...
            "ngsolve_capability_checked", true);

    case "galerkin_ngsolve"
        fixture = options.Fixture;
        if fixture == ""
            fixture = "unit_sphere_coarse.vol";
        end
        volPath = fullfile(root, "fixtures", "mesh_topology", fixture);
        raw = verifyGalerkinAgainstNgsolve(volPath);
        report.reference = "ngsolve_bem_laplace_galerkin_reference";
        report.fixture = fixture;
        report.ngsolve_version = raw.ngsolveVersion;
        report.reference_intorder = raw.referenceIntorder;
        report.mass_rel_diff = raw.massRelDiff;
        report.operator_v_rel_diff = raw.operatorVRelDiff;
        report.operator_k_rel_diff = raw.operatorKRelDiff;
        report.capacitance_rel_diff = raw.capacitanceRelDiff;
        report.relative_error = max([raw.operatorVRelDiff, ...
            raw.operatorKRelDiff, raw.capacitanceRelDiff]);
        report.checks = raw.checks;
        report.checks.max_error_within_tolerance = report.relative_error <= tol;

    case "helmholtz_ngsolve"
        fixture = options.Fixture;
        if fixture == ""
            fixture = "unit_sphere_coarse.vol";
        end
        volPath = fullfile(root, "fixtures", "mesh_topology", fixture);
        raw = verifyHelmholtzAgainstNgsolve(volPath);
        report.reference = "ngsolve_bem_helmholtz_three_way_reference";
        report.fixture = fixture;
        report.ngsolve_version = raw.ngsolveVersion;
        errors = collectHelmholtzErrors(raw);
        report.relative_error = max(errors);
        report.num_cases = numel(raw.cases);
        report.case_errors = errors;
        report.checks = struct( ...
            "all_cases_ok", all(arrayfun(@(x) x.status == "ok", raw.cases)), ...
            "max_error_within_tolerance", report.relative_error <= tol);
end

report.duration_s = toc(timer);
report.pass = all(structfun(@(x) logical(x), report.checks)) ...
    && isfinite(report.relative_error) ...
    && report.relative_error <= tol;
if report.pass
    report.status = "ok";
end
end


function errors = collectHelmholtzErrors(raw)
errors = zeros(0, 1);
for k = 1:numel(raw.cases)
    c = raw.cases(k);
    errors(end + 1, 1) = c.operatorRelDiff; %#ok<AGROW>
    errors(end + 1, 1) = c.doubleLayerRelDiff; %#ok<AGROW>
    errors(end + 1, 1) = c.pointSourceDensityRelDiff; %#ok<AGROW>
    errors(end + 1, 1) = c.planeWaveDensityRelDiff; %#ok<AGROW>
    errors(end + 1, 1) = c.pointSourceProbeCrossCode; %#ok<AGROW>
    errors(end + 1, 1) = c.planeWaveProbeCrossCode; %#ok<AGROW>
    errors(end + 1, 1) = c.rigidTraceCrossCode; %#ok<AGROW>
    errors(end + 1, 1) = c.rigidProbeCrossCode; %#ok<AGROW>
end
end
