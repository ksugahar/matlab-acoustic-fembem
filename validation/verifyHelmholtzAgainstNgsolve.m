function report = verifyHelmholtzAgainstNgsolve(volFile, matFile)
%VERIFYHELMHOLTZAGAINSTNGSOLVE 3-way Helmholtz check: analytic / here / ngbem.
%
%   report = verifyHelmholtzAgainstNgsolve( ...
%       "fixtures/mesh_topology/unit_sphere_fine.vol");
%
% The acoustic (Helmholtz single-layer) leg of the cross-validation ladder,
% run as a THREE-way comparison on the same .vol mesh and the same
% continuous-P1 space:
%
%   analytic   interior point source (exact, no truncation) and the
%              sound-soft sphere partial-wave series
%   this repo  GalerkinSingleLayer Wavenumber path + singleLayerDirichletSolve
%   ngsolve    ngsolve.bem HelmholtzSingleLayerPotentialOperator reference
%              (committed .mat from exportNgsolveBemHelmholtzReference.py,
%              with NGSolve's own GetPotential probe values)
%
% The intended reading, locked from the 2026-07-03 measurements: the two
% CODES agree with each other 10-30x tighter than either agrees with the
% analytic SPHERE (probe cross-code 3e-4..6e-3 vs 1-10% against analytic),
% which proves the analytic deviation is the faceted-geometry
% discretization error, not an implementation bug. The kernel convention is
% pinned per case: our operator must match ngbem far better than its
% conjugate (e^{+ikr} on both sides).

arguments
    volFile (1,1) string
    matFile (1,1) string = ""
end

if matFile == ""
    [~, base] = fileparts(volFile);
    matFile = fullfile(fileparts(mfilename("fullpath")), "data", ...
        "ngbem_helmholtz_reference_" + base + ".mat");
end
if ~isfile(matFile)
    error("verifyHelmholtzAgainstNgsolve:reference", ...
        "ngsolve.bem Helmholtz reference not found: %s\nRegenerate with:\n  python validation/exportNgsolveBemHelmholtzReference.py %s %s", ...
        matFile, volFile, matFile);
end
S = load(matFile);

mesh = VolMesh(volFile);
surface = mesh.boundary();
ids = surface.volNodeIds;

vertexMismatch = max(abs(surface.vtx - S.vtx(ids, :)), [], "all");
if vertexMismatch > 1e-12
    error("verifyHelmholtzAgainstNgsolve:vertices", ...
        "NGSolve vertex order does not match the .vol point order (max %.3e).", ...
        vertexMismatch);
end

sourcePoint = S.sourcePoint(:).';
probes = S.probePoints;
caseNames = string(fieldnames(S));
caseNames = caseNames(startsWith(caseNames, "case_"));

report = struct();
report.kind = "helmholtz_single_layer_three_way_cross_check";
report.volFile = string(volFile);
report.matFile = string(matFile);
report.nNodes = size(surface.vtx, 1);
report.ngsolveVersion = string(S.ngsolveVersion);
report.quadratureOrder = 7;
report.cases = struct([]);

for c = 1:numel(caseNames)
    C = S.(caseNames(c));
    k = C.k;
    r = struct();
    r.name = caseNames(c);
    r.wavenumber = k;
    r.referenceIntorderConvergenceV = C.intorderConvergenceV;

    Vng = C.V(ids, ids);
    Vours = GalerkinSingleLayer(surface, "Wavenumber", k, ...
        "QuadratureOrder", 7).matrix;
    r.operatorRelDiff = norm(Vours - Vng, "fro") / norm(Vng, "fro");
    r.operatorConjRelDiff = norm(Vours - conj(Vng), "fro") / norm(Vng, "fro");

    % Helmholtz double layer leg (measured: 1.8-3.1e-3, conj 0.12-1.4)
    Kng = C.K(ids, ids);
    Kours = GalerkinDoubleLayer(surface, "Wavenumber", k, ...
        "QuadratureOrder", 7).matrix;
    r.doubleLayerRelDiff = norm(Kours - Kng, "fro") / norm(Kng, "fro");
    r.doubleLayerConjRelDiff = norm(Kours - conj(Kng), "fro") / norm(Kng, "fro");
    r.referenceIntorderConvergenceK = C.intorderConvergenceK;

    % rigid (sound-hard) leg: total-field K equation, second-kind - the
    % cross-code agreement (1e-5..1e-4) is 1-2 orders TIGHTER than the
    % first-kind V solves, the second-kind conditioning made visible
    rig = rigidScatteringSolve(surface, "Wavenumber", k, ...
        "QuadratureOrder", 7);
    tng = C.tRigid(ids);
    r.rigidTraceCrossCode = norm(rig.trace - tng(:)) / norm(tng(:));
    refRig = rigidSphereScattering(k, 1.0, probes);
    usRig = rig.scatteredAt(probes);
    pngRig = C.probeRigidScattered(:);
    r.rigidProbeCrossCode = max(abs(usRig - pngRig) ./ abs(pngRig));
    r.rigidProbeVsAnalyticOurs = ...
        max(abs(usRig - refRig.scattered) ./ abs(refRig.scattered));
    r.rigidProbeVsAnalyticNgsolve = ...
        max(abs(pngRig - refRig.scattered) ./ abs(refRig.scattered));

    % point-source leg (exact analytic gate)
    g = acousticPointSource(k, sourcePoint, surface.vtx);
    gng = C.gPointSource(ids);
    r.pointSourceDataIdentity = max(abs(g - gng(:)));
    sol = singleLayerDirichletSolve(surface, g, ...
        "Wavenumber", k, "QuadratureOrder", 7);
    qng = C.qPointSource(ids);
    r.pointSourceDensityRelDiff = norm(sol.q - qng(:)) / norm(qng(:));
    pours = sol.potentialAt(probes);
    png = C.probePointSource(:);
    pana = acousticPointSource(k, sourcePoint, probes);
    r.pointSourceProbeCrossCode = max(abs(pours - png) ./ abs(png));
    r.pointSourceProbeVsAnalyticOurs = max(abs(pours - pana) ./ abs(pana));
    r.pointSourceProbeVsAnalyticNgsolve = max(abs(png - pana) ./ abs(pana));

    % plane-wave (sound-soft scattering) leg vs the partial-wave series
    gs = -exp(1i * k * surface.vtx(:, 3));
    gsng = C.gPlaneWave(ids);
    r.planeWaveDataIdentity = max(abs(gs - gsng(:)));
    sols = singleLayerDirichletSolve(surface, gs, ...
        "Wavenumber", k, "QuadratureOrder", 7);
    qsng = C.qPlaneWave(ids);
    r.planeWaveDensityRelDiff = norm(sols.q - qsng(:)) / norm(qsng(:));
    ps = sols.potentialAt(probes);
    psng = C.probePlaneWave(:);
    ref = softSphereScattering(k, 1.0, probes);
    r.seriesTruncationTail = ref.truncationTail;
    r.planeWaveProbeCrossCode = max(abs(ps - psng) ./ abs(psng));
    r.planeWaveProbeVsSeriesOurs = max(abs(ps - ref.scattered) ./ abs(ref.scattered));
    r.planeWaveProbeVsSeriesNgsolve = max(abs(psng - ref.scattered) ./ abs(ref.scattered));

    r.checks = struct( ...
        "referenceConverged", ...
            C.intorderConvergenceV < 1e-6 && C.intorderConvergenceK < 1e-6, ...
        "boundaryDataIdentical", ...
            r.pointSourceDataIdentity < 1e-12 && r.planeWaveDataIdentity < 1e-12, ...
        "conventionPinned", ...
            r.operatorConjRelDiff > 0.1 && r.operatorRelDiff < 0.1 * r.operatorConjRelDiff, ...
        "operatorClose", r.operatorRelDiff < 2e-2, ...
        "doubleLayerClose", r.doubleLayerRelDiff < 1e-2, ...
        "doubleLayerConventionPinned", ...
            r.doubleLayerConjRelDiff > 0.05 && ...
            r.doubleLayerRelDiff < 0.1 * r.doubleLayerConjRelDiff, ...
        "densityClose", ...
            r.pointSourceDensityRelDiff < 2e-2 && r.planeWaveDensityRelDiff < 2e-2, ...
        "probeCrossCodeClose", ...
            r.pointSourceProbeCrossCode < 2e-2 && r.planeWaveProbeCrossCode < 2e-2, ...
        "rigidCrossCodeClose", ...
            r.rigidTraceCrossCode < 1e-3 && r.rigidProbeCrossCode < 1e-3, ...
        "probeAnalyticWithinFaceting", ...
            r.pointSourceProbeVsAnalyticOurs < 0.15 && ...
            r.planeWaveProbeVsSeriesOurs < 0.15 && ...
            r.pointSourceProbeVsAnalyticNgsolve < 0.15 && ...
            r.planeWaveProbeVsSeriesNgsolve < 0.15, ...
        "rigidWithinFaceting", ...
            r.rigidProbeVsAnalyticOurs < 0.20 && ...
            r.rigidProbeVsAnalyticNgsolve < 0.20, ...
        "seriesTailNegligible", r.seriesTruncationTail < 1e-12);
    r.status = "needs_attention";
    if all(structfun(@(x) logical(x), r.checks))
        r.status = "ok";
    end
    if isempty(report.cases)
        report.cases = r;
    else
        report.cases(end + 1) = r; %#ok<AGROW>
    end
end

if all(arrayfun(@(x) x.status == "ok", report.cases))
    report.status = "ok";
else
    report.status = "needs_attention";
end
end
