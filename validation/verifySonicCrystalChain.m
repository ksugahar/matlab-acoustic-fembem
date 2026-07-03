function report = verifySonicCrystalChain(volFile, matFile)
%VERIFYSONICCRYSTALCHAIN 4-leg check of the multi-body sphere-chain rung.
%
%   report = verifySonicCrystalChain();
%
% The sonic-crystal chain fixture (five sound-soft spheres, R = 0.3,
% d = 1.5) checked FOUR ways at k = 1.0 (sub-wavelength regime) and
% k = 2.0944 (the Bragg wavenumber pi/d of the chain):
%
%   exact analytic   interior point source inside the middle sphere ->
%                    exterior reproduction (multi-body uniqueness gate)
%   analytic class   Foldy monopole multiple scattering (low-k reference;
%                    its l >= 1 truncation grows with k and is REPORTED,
%                    not hidden)
%   ngsolve.bem      same-mesh continuous-P1 Helmholtz V + GetPotential
%                    probe values (committed .mat)
%   NGSolve FEM      volume Helmholtz solve on an entirely different
%                    discretization (air ball, Dirichlet spheres, order-2,
%                    first-order Sommerfeld ABC) at the same probes
%
% Physics finding locked by the tests (see testSonicCrystalChain): the
% sparse free-space chain shows broadband sub-wavelength attenuation and
% NO Bragg stop band - band gaps need confinement or transverse
% periodicity (the duct/Bloch FEM rung).
%
% MATLAB operators are assembled at gss 3 here (the chain has 634
% triangles; the 7-point smooth-correction kernel would be a ~315 MB
% transient) - the bands below are locked for gss 3.

arguments
    volFile (1,1) string = defaultVolFile()
    matFile (1,1) string = defaultMatFile()
end

if ~isfile(matFile)
    error("verifySonicCrystalChain:reference", ...
        "NGSolve chain reference not found: %s\nRegenerate with:\n  python validation/exportNgsolveChainReference.py", ...
        matFile);
end
S = load(matFile);

mesh = VolMesh(volFile);
surface = mesh.boundary();
ids = surface.volNodeIds;

vertexMismatch = max(abs(surface.vtx - S.vtx(ids, :)), [], "all");
if vertexMismatch > 1e-12
    error("verifySonicCrystalChain:vertices", ...
        "NGSolve vertex order does not match the .vol point order (max %.3e).", ...
        vertexMismatch);
end

sourcePoint = S.sourcePoint(:).';
probes = S.probePoints;
radius = S.radius;
centers = S.centers;
caseNames = string(fieldnames(S));
caseNames = caseNames(startsWith(caseNames, "case_"));

report = struct();
report.kind = "sonic_crystal_chain_four_leg_cross_check";
report.volFile = string(volFile);
report.matFile = string(matFile);
report.nNodes = size(surface.vtx, 1);
report.ngsolveVersion = string(S.ngsolveVersion);
report.quadratureOrder = 3;
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
        "QuadratureOrder", 3).matrix;
    r.operatorRelDiff = norm(Vours - Vng, "fro") / norm(Vng, "fro");
    r.operatorConjRelDiff = norm(Vours - conj(Vng), "fro") / norm(Vng, "fro");

    % exact multi-body gate: point source inside the middle sphere
    g = acousticPointSource(k, sourcePoint, surface.vtx);
    sol = singleLayerDirichletSolve(surface, g, ...
        "Wavenumber", k, "QuadratureOrder", 3);
    qng = C.qPointSource(ids);
    r.pointSourceDensityRelDiff = norm(sol.q - qng(:)) / norm(qng(:));
    pours = sol.potentialAt(probes);
    pana = acousticPointSource(k, sourcePoint, probes);
    png = C.probePointSource(:);
    r.pointSourceExactGate = max(abs(pours - pana) ./ abs(pana));
    r.pointSourceProbeCrossCode = max(abs(pours - png) ./ abs(png));

    % plane-wave leg: BEM (ours) vs ngbem vs FEM vs Foldy
    gs = -exp(1i * k * surface.vtx(:, 3));
    sols = singleLayerDirichletSolve(surface, gs, ...
        "Wavenumber", k, "QuadratureOrder", 3);
    qsng = C.qPlaneWave(ids);
    r.planeWaveDensityRelDiff = norm(sols.q - qsng(:)) / norm(qsng(:));
    pinc = exp(1i * k * probes(:, 3));
    ptot = sols.potentialAt(probes) + pinc;
    ptotNg = C.probePlaneWave(:) + pinc;
    ptotFem = C.femProbePlaneWave(:);
    foldy = foldyPointScattering(k, radius, centers, probes);
    r.planeWaveProbeCrossCode = max(abs(ptot - ptotNg) ./ abs(ptotNg));
    r.planeWaveProbeVsFem = max(abs(ptot - ptotFem) ./ abs(ptotFem));
    r.planeWaveProbeVsFoldy = max(abs(ptot - foldy.total) ./ abs(foldy.total));

    r.checks = struct( ...
        "referenceConverged", C.intorderConvergenceV < 1e-6, ...
        "conventionPinned", ...
            r.operatorConjRelDiff > 0.1 && ...
            r.operatorRelDiff < 0.1 * r.operatorConjRelDiff, ...
        "operatorClose", r.operatorRelDiff < 3e-2, ...
        "pointSourceExact", r.pointSourceExactGate < 4e-2, ...
        "probeCrossCodeClose", ...
            r.pointSourceProbeCrossCode < 3e-2 && ...
            r.planeWaveProbeCrossCode < 3e-2, ...
        "femLegAgrees", r.planeWaveProbeVsFem < 0.10, ...
        "foldyReferenceBehaves", r.planeWaveProbeVsFoldy < foldyBand(k));
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


function band = foldyBand(k)
%FOLDYBAND Honest per-k band for the monopole-only Foldy reference.
% Measured single-sphere l >= 1 truncation: ~9% at kR = 0.3, ~31% at
% kR = 0.66; the chain probes sit in the forward direction where the
% collective low-k agreement is tighter (~3% at k = 1.0).
if k <= 1.2
    band = 0.08;
else
    band = 0.35;
end
end


function file = defaultVolFile()
file = fullfile(fileparts(mfilename("fullpath")), "..", ...
    "fixtures", "mesh_topology", "soft_sphere_chain_5.vol");
end


function file = defaultMatFile()
file = fullfile(fileparts(mfilename("fullpath")), "data", ...
    "ngsolve_chain_reference_soft_sphere_chain_5.mat");
end
