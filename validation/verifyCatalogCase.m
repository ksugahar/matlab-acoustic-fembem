function result = verifyCatalogCase(item, caps)
%VERIFYCATALOGCASE Verify one non-mesh catalog case with readable gates.
%
% The gates are intentionally small. They make sure the MATLAB teaching
% operators agree with NGSolve where a direct first-order reference is
% available, and otherwise check the analytic identities that students should
% be able to inspect in the source.

arguments
    item (1,1) struct
    caps (1,1) struct
end

failures = strings(0, 1);
details = struct();

try
    failures = requirePass(failures, caps.ok, "NGSolve capability query did not pass.");
    switch string(item.category)
        case "02_h1_scalar_fem"
            [failures, details] = verifyH1Scalar(item, caps, failures);
        case "03_hcurl_edge_fem"
            [failures, details] = verifyHCurlEdge(item, caps, failures);
        case "04_laplace_bem_dense"
            [failures, details] = verifyLaplaceDense(item, caps, failures);
        case "05_laplace_hmatrix"
            [failures, details] = verifyLaplaceHMatrix(item, caps, failures);
        case "06_acoustic_low_frequency"
            [failures, details] = verifyAcousticLowFrequency(item, caps, failures);
        case "07_acoustic_helmholtz"
            [failures, details] = verifyAcousticHelmholtz(item, caps, failures);
        case "08_scalar_fem_bem_coupling"
            [failures, details] = verifyScalarCoupling(item, caps, failures);
        case "09_rwg_hcurl_trace"
            [failures, details] = verifyRwgHCurlTrace(item, caps, failures);
        case "10_ngsolve_bem_reference"
            [failures, details] = verifyNgsolveBemReference(item, caps, failures);
        otherwise
            failures(end + 1, 1) = "No verifier for category " + string(item.category);
    end
catch err
    failures(end + 1, 1) = string(err.identifier) + ": " + string(err.message);
end

result = struct();
result.id = item.id;
result.title = item.title;
result.category = item.category;
result.reference = item.reference;
result.tolerance = item.tolerance;
result.passed = isempty(failures);
result.failures = failures;
result.details = details;
end


function [failures, details] = verifyH1Scalar(item, caps, failures)
model = unitModel();
[~, stiff] = model.h1.stiffness();
[~, mass] = model.h1.mass();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "h1_scalar_fem", "slot", slot);
switch slot
    case 1
        err = norm(full(stiff.stiffness) - caps.h1Stiffness, "fro");
        details.error = err;
        failures = requirePass(failures, err < tol, "P1 stiffness differs from NGSolve.");
    case 2
        err = norm(full(mass.mass) - caps.h1Mass, "fro");
        details.error = err;
        failures = requirePass(failures, err < tol, "P1 mass differs from NGSolve.");
    case 3
        rowSum = norm(full(stiff.stiffness) * ones(4, 1));
        details.rowSum = rowSum;
        failures = requirePass(failures, rowSum < tol, "Constant field should have zero stiffness residual.");
    case 4
        u = sum(model.mesh.vtx, 2);
        energy = u.' * stiff.stiffness * u;
        details.energy = full(energy);
        failures = requirePass(failures, abs(energy - 0.5) < tol, ...
            "Linear manufactured potential energy should equal volume*|grad u|^2.");
    case 5
        Kdd = full(stiff.stiffness(2:4, 2:4));
        lambda = eig(Kdd);
        details.minEigenvalue = min(lambda);
        failures = requirePass(failures, all(lambda > 0), "Dirichlet eliminated stiffness is not positive definite.");
    case 6
        face = [1, 2, 3];
        faceArea = triangleArea(model.mesh.vtx(face, :));
        load = zeros(4, 1);
        load(face) = faceArea / 3;
        details.faceLoadSum = sum(load);
        failures = requirePass(failures, abs(sum(load) - 0.5) < tol, "Neumann face load sum is wrong.");
    case 7
        patch = fourTetModel();
        [~, K] = patch.h1.stiffness();
        details.nodes = size(patch.mesh.vtx, 1);
        details.tets = size(patch.mesh.tet, 1);
        failures = requirePass(failures, issymmetricWithin(K.stiffness, tol), "Two-tet patch stiffness is not symmetric.");
        failures = requirePass(failures, norm(full(K.stiffness) * ones(size(K.stiffness, 1), 1)) < tol, ...
            "Two-tet patch does not preserve constants.");
    case 8
        [~, scaled] = model.h1.stiffness(2.0);
        err = norm(full(scaled.stiffness - 2.0 * stiff.stiffness), "fro");
        details.error = err;
        failures = requirePass(failures, err < tol, "Material coefficient scaling failed for H1 stiffness.");
    case 9
        residual = full(stiff.stiffness) * ones(4, 1);
        details.residualNorm = norm(residual);
        failures = requirePass(failures, norm(residual) < tol, "Volume residual balance failed for constants.");
    case 10
        errK = norm(full(stiff.stiffness) - caps.h1Stiffness, "fro");
        errM = norm(full(mass.mass) - caps.h1Mass, "fro");
        details.stiffnessError = errK;
        details.massError = errM;
        failures = requirePass(failures, caps.h1Dofs == 4, "NGSolve H1 dof count changed.");
        failures = requirePass(failures, errK < tol && errM < tol, "Combined NGSolve H1 comparison failed.");
end
end


function [failures, details] = verifyHCurlEdge(item, caps, failures)
model = unitModel();
[~, ~, ned] = model.hcurl.matrices();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "hcurl_edge_fem", "slot", slot);
switch slot
    case 1
        details.edges = size(ned.edges, 1);
        failures = requirePass(failures, size(ned.edges, 1) == 6 && caps.hcurlDofs == 6, ...
            "Single tetrahedron should have six Nedelec0 edge dofs.");
    case 2
        details.signs = model.hcurl.tetEdgeSigns;
        failures = requirePass(failures, all(abs(model.hcurl.tetEdgeSigns(:)) == 1), ...
            "HCurl edge signs must be +/-1.");
    case 3
        patch = fourTetModel();
        nTetEdges = 6 * size(patch.mesh.tet, 1);
        details.uniqueEdges = size(patch.hcurl.edges, 1);
        details.rawTetEdges = nTetEdges;
        failures = requirePass(failures, size(patch.hcurl.edges, 1) < nTetEdges, ...
            "Shared-edge continuity did not reduce global HCurl dofs.");
    case 4
        curlNorms = sqrt(sum(ned.localCurls(:, :, 1).^2, 2));
        details.minCurlNorm = min(curlNorms);
        failures = requirePass(failures, all(curlNorms > 0), "Nedelec0 basis curls should be nonzero on a tetrahedron.");
    case 5
        err = norm(full(ned.mass - ned.mass.'), "fro");
        details.symmetryError = err;
        failures = requirePass(failures, err < tol, "Nedelec0 mass matrix is not symmetric.");
    case 6
        err = norm(full(ned.curlCurl - ned.curlCurl.'), "fro");
        details.symmetryError = err;
        failures = requirePass(failures, err < tol, "Nedelec0 curl-curl matrix is not symmetric.");
    case 7
        G = edgeGradientMatrix(ned.edges, size(model.mesh.vtx, 1));
        err = norm(full(ned.curlCurl * G), "fro");
        details.gradientCurlError = err;
        failures = requirePass(failures, err < tol, "Gradient edge fields should be in the curl-curl nullspace.");
    case 8
        details.boundaryEdges = numel(model.rwg.dofEdgeIds);
        failures = requirePass(failures, numel(model.rwg.dofEdgeIds) == 6, ...
            "Closed tetrahedron boundary should expose six RWG edge dofs.");
    case 9
        [~, ~, scaled] = model.hcurl.matrices(3.0);
        err = norm(full(scaled.mass - 3.0 * ned.mass), "fro") + ...
            norm(full(scaled.curlCurl - 3.0 * ned.curlCurl), "fro");
        details.error = err;
        failures = requirePass(failures, err < tol, "Material scaling failed for HCurl matrices.");
    case 10
        errMass = norm(full(ned.mass) - caps.hcurlMass, "fro");
        errCurl = norm(full(ned.curlCurl) - caps.hcurlCurlCurl, "fro");
        details.massError = errMass;
        details.curlCurlError = errCurl;
        failures = requirePass(failures, errMass < tol && errCurl < tol, "NGSolve HCurl matrix comparison failed.");
end
end


function [failures, details] = verifyLaplaceDense(item, caps, failures)
model = unitModel();
[centers, weights] = triangleCentroidsAndAreas(model.surface.vtx, model.surface.tri);
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "laplace_bem_dense", "slot", slot);
switch slot
    case 1
        target = [0, 0, 0];
        source = [1, 0, 0];
        K = HelmholtzKernel(target, source);
        details.value = K.singleLayer(1, 1);
        failures = requirePass(failures, abs(K.singleLayer(1, 1) - 1 / (4 * pi)) < tol, ...
            "Single-point Laplace kernel is wrong.");
    case 2
        source = centers + [2, 0, 0];
        K = HelmholtzKernel(centers, source, "SourceWeights", weights);
        details.matrixSize = size(K.singleLayer);
        failures = requirePass(failures, all(K.singleLayer(:) > 0), "Separated Laplace kernel should be positive.");
    case 3
        K = HelmholtzKernel(centers, centers, "SourceWeights", weights);
        weighted = diag(weights) * K.singleLayer;
        err = norm(weighted - weighted.', "fro");
        details.weightedSymmetryError = err;
        failures = requirePass(failures, err < tol, "Weighted Laplace single-layer matrix is not symmetric.");
    case 4
        K = HelmholtzKernel(centers, centers, "DiagonalValue", 0.125);
        details.diagonal = diag(K.singleLayer).';
        failures = requirePass(failures, max(abs(diag(K.singleLayer) - 0.125)) < tol, ...
            "Near-field diagonal policy was not applied.");
    case 5
        sphere = spherePoints(12);
        K = HelmholtzKernel(sphere, sphere + [0.25, 0, 0], "DiagonalValue", 0.0);
        details.meanPotential = mean(real(K.singleLayer * ones(size(sphere, 1), 1)));
        failures = requirePass(failures, details.meanPotential > 0, "Coarse sphere Laplace potential should be positive.");
    case 6
        cube = cubeSurfacePoints();
        K = HelmholtzKernel(cube, cube + [0.25, 0, 0], "DiagonalValue", 0.0);
        details.meanPotential = mean(real(K.singleLayer * ones(size(cube, 1), 1)));
        failures = requirePass(failures, details.meanPotential > 0, "Coarse cube Laplace potential should be positive.");
    case 7
        K = HelmholtzKernel(centers, centers, "SourceWeights", weights, "DiagonalValue", 0.0);
        potential = K.singleLayer * ones(size(centers, 1), 1);
        details.minPotential = min(potential);
        failures = requirePass(failures, all(potential >= 0), "Constant density Laplace potential should be nonnegative.");
    case 8
        K = HelmholtzKernel([0, 0, 1], [0, 0, 0], "SourceNormals", [0, 0, 1]);
        details.doubleLayerValue = K.doubleLayerSourceNormal(1, 1);
        failures = requirePass(failures, K.doubleLayerSourceNormal(1, 1) > 0, ...
            "Source-normal double-layer sign is unexpected.");
    case 9
        K = HelmholtzKernel(centers, centers + [2, 0, 0], "SourceWeights", weights);
        x = (1:size(centers, 1)).';
        y1 = K.singleLayer * x;
        y2 = K.singleLayer * x;
        details.reproducibilityError = norm(y1 - y2);
        failures = requirePass(failures, norm(y1 - y2) < tol, "Dense Laplace matvec is not reproducible.");
    case 10
        failures = requirePass(failures, caps.hasBem && caps.hasLaplaceSL, "NGSolve Laplace BEM is not available.");
        K0 = HelmholtzKernel(centers, centers + [2, 0, 0]);
        details.referenceNorm = norm(K0.singleLayer, "fro");
        failures = requirePass(failures, isfinite(details.referenceNorm) && details.referenceNorm > 0, ...
            "Laplace BEM reference matrix norm is invalid.");
end
end


function [failures, details] = verifyLaplaceHMatrix(item, caps, failures)
[target, source] = separatedPointClouds();
near = target;
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "laplace_hmatrix", "slot", slot);
switch slot
    case 1
        H = HMatrix(near, [], "LeafSize", 2, "Eta", 1.0, "RankTolerance", 1e-12);
        stats = H.stats();
        details.splitBlocks = stats.splitBlocks;
        failures = requirePass(failures, stats.splitBlocks > 0, "Cluster tree did not split.");
    case 2
        Hsmall = HMatrix(target, source, "LeafSize", 2, "Eta", 0.2);
        Hlarge = HMatrix(target, source, "LeafSize", 2, "Eta", 2.0);
        sSmall = Hsmall.stats();
        sLarge = Hlarge.stats();
        details.blocksSmallEta = sSmall.blocks;
        details.blocksLargeEta = sLarge.blocks;
        details.lowRankSmallEta = sSmall.lowRankBlocks;
        details.lowRankLargeEta = sLarge.lowRankBlocks;
        failures = requirePass(failures, sLarge.blocks <= sSmall.blocks && sLarge.lowRankBlocks > 0, ...
            "Larger admissibility eta should accept coarser far blocks.");
    case 3
        H = HMatrix(target, source, "LeafSize", 8, "Eta", 2.0);
        stats = H.stats();
        details.lowRankBlocks = stats.lowRankBlocks;
        failures = requirePass(failures, stats.lowRankBlocks > 0, "Separated far field was not compressed.");
    case 4
        H = HMatrix(near, [], "LeafSize", 64, "Eta", 0.1);
        stats = H.stats();
        details.denseBlocks = stats.denseBlocks;
        failures = requirePass(failures, stats.denseBlocks == 1, "Near block should stay dense.");
    case 5
        H = HMatrix(target, source, "LeafSize", 2, "Eta", 2.0, "RankTolerance", 1e-14);
        x = sin((1:size(source, 1)).');
        yH = H.matvec(x);
        K = HelmholtzKernel(target, source);
        err = norm(yH - K.singleLayer * x);
        details.matvecError = err;
        failures = requirePass(failures, err < tol, "H-matrix matvec differs from dense Laplace matvec.");
    case 6
        Hlo = HMatrix(target, source, "RankTolerance", 1e-12);
        Hhi = HMatrix(target, source, "RankTolerance", 1e-2);
        slo = Hlo.stats();
        shi = Hhi.stats();
        details.maxRankLowTolerance = slo.maxRank;
        details.maxRankHighTolerance = shi.maxRank;
        failures = requirePass(failures, slo.maxRank >= shi.maxRank, "Rank tolerance sweep is not monotone.");
    case 7
        Hsmall = HMatrix(near, [], "LeafSize", 2);
        Hlarge = HMatrix(near, [], "LeafSize", 8);
        sSmall = Hsmall.stats();
        sLarge = Hlarge.stats();
        details.blocksSmallLeaf = sSmall.blocks;
        details.blocksLargeLeaf = sLarge.blocks;
        failures = requirePass(failures, sSmall.blocks >= sLarge.blocks, "Smaller leaf size should create more blocks.");
    case 8
        H = HMatrix(target, source, "LeafSize", 2);
        stats = H.stats();
        details.compressionRatio = stats.compressionRatio;
        failures = requirePass(failures, isfinite(stats.compressionRatio) && stats.compressionRatio > 0, ...
            "H-matrix storage ratio is invalid.");
    case 9
        sphere = spherePoints(16);
        H = HMatrix(sphere, [], "LeafSize", 2, "DiagonalValue", 0.0);
        y = H.matvec(ones(size(sphere, 1), 1));
        details.meanPotential = mean(y);
        failures = requirePass(failures, all(isfinite(y)) && mean(y) > 0, "Sphere H-matrix Laplace response is invalid.");
    case 10
        H = HMatrix(target, source, "LeafSize", 2);
        stats = H.stats();
        details.blocks = stats.blocks;
        failures = requirePass(failures, caps.hasBem && caps.hasLaplaceSL && stats.blocks > 0, ...
            "NGSolve BEM/H-matrix comparison gate is unavailable.");
end
end


function [failures, details] = verifyAcousticLowFrequency(item, caps, failures)
[target, source] = separatedPointClouds();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "acoustic_low_frequency", "slot", slot);
switch slot
    case 1
        K0 = HelmholtzKernel(target, source, "Wavenumber", 0);
        L = HelmholtzKernel(target, source);
        err = norm(K0.singleLayer - L.singleLayer, "fro");
        details.error = err;
        failures = requirePass(failures, err < tol, "Helmholtz k=0 limit does not match Laplace.");
    case 2
        k = 1e-7;
        K = HelmholtzKernel([0, 0, 0], [1, 0, 0], "Wavenumber", k);
        expected = 1 / (4 * pi) + 1i * k / (4 * pi);
        details.error = abs(K.singleLayer(1, 1) - expected);
        failures = requirePass(failures, details.error < 1e-12, "expm1 single-layer correction is inaccurate.");
    case 3
        k = 1e-8;
        K = HelmholtzKernel(target, source, "Wavenumber", k, "TaylorCutoff", 1);
        direct = exp(1i * k * pairwiseDistance(target, source)) ./ (4 * pi * pairwiseDistance(target, source));
        err = norm(K.singleLayer - direct, "fro");
        details.error = err;
        failures = requirePass(failures, err < 1e-12, "Taylor single-layer correction is inaccurate.");
    case 4
        K = HelmholtzKernel([0, 0, 1], [0, 0, 0], ...
            "SourceNormals", [0, 0, 1], "Wavenumber", 1e-8, "TaylorCutoff", 1);
        details.correction = K.doubleLayerSourceNormalCorrection(1, 1);
        failures = requirePass(failures, abs(K.doubleLayerSourceNormalCorrection(1, 1)) < 1e-14, ...
            "Low-frequency double-layer correction is not small.");
    case 5
        ks = [0, 1e-8, 1e-6, 1e-4];
        norms = zeros(size(ks));
        for k = 1:numel(ks)
            K = HelmholtzKernel(target, source, "Wavenumber", ks(k));
            norms(k) = norm(K.singleLayer, "fro");
        end
        details.norms = norms;
        failures = requirePass(failures, all(isfinite(norms)) && max(abs(diff(norms))) < 1e-3, ...
            "Small-kr sweep is not smooth.");
    case 6
        pts = target;
        K = HelmholtzKernel(pts, pts, "Wavenumber", 1e-5, "DiagonalValue", 0);
        err = norm(K.singleLayer - K.singleLayer.', "fro");
        details.complexSymmetryError = err;
        failures = requirePass(failures, err < tol, "Low-frequency single-layer matrix is not complex symmetric.");
    case 7
        K = HelmholtzKernel(target, source, "Wavenumber", 1e-6);
        x = cos((1:size(source, 1)).');
        y = K.singleLayer * x;
        details.responseNorm = norm(y);
        failures = requirePass(failures, all(isfinite(y)) && norm(y) > 0, "Low-frequency matvec is invalid.");
    case 8
        sphere = spherePoints(12);
        K = HelmholtzKernel(sphere, sphere + [2, 0, 0], "Wavenumber", 1e-6);
        details.monopoleMean = mean(real(K.singleLayer * ones(size(sphere, 1), 1)));
        failures = requirePass(failures, details.monopoleMean > 0, "Low-frequency sphere monopole response is invalid.");
    case 9
        H = HMatrix(target, source, "RankTolerance", 1e-12);
        K = HelmholtzKernel(target, source, "Wavenumber", 1e-8);
        x = ones(size(source, 1), 1);
        err = norm(H.matvec(x) - real(K.singleLayer * x));
        details.hmatrixCandidateError = err;
        failures = requirePass(failures, err < 1e-7, "Low-frequency H-matrix candidate differs from Laplace limit.");
    case 10
        K = HelmholtzKernel(target, source, "Wavenumber", 1e-8);
        details.referenceNorm = norm(K.singleLayer, "fro");
        failures = requirePass(failures, caps.hasBem && caps.hasHelmholtzSL && isfinite(details.referenceNorm), ...
            "NGSolve low-frequency BEM comparison gate is unavailable.");
end
end


function [failures, details] = verifyAcousticHelmholtz(item, caps, failures)
[target, source] = separatedPointClouds();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "acoustic_helmholtz", "slot", slot);
switch slot
    case 1
        K = HelmholtzKernel([0, 0, 0], [2, 0, 0], "Wavenumber", 3);
        expected = exp(1i * 6) / (8 * pi);
        details.error = abs(K.singleLayer(1, 1) - expected);
        failures = requirePass(failures, details.error < tol, "Two-point Helmholtz kernel is wrong.");
    case 2
        k = 4;
        plane = exp(1i * k * target(:, 1));
        details.maxAmplitudeError = max(abs(abs(plane) - 1));
        failures = requirePass(failures, details.maxAmplitudeError < tol, "Plane-wave boundary data lost unit amplitude.");
    case 3
        sphere = spherePoints(12);
        op = AcousticSingleLayer(sphere, sphere + [2, 0, 0], "Wavenumber", 0.1);
        details.responseNorm = norm(op.apply(ones(size(op.matrix, 2), 1)));
        failures = requirePass(failures, isfinite(details.responseNorm) && details.responseNorm > 0, ...
            "Low-ka pulsating sphere toy response is invalid.");
    case 4
        sphere = spherePoints(12);
        op = AcousticSingleLayer(sphere, sphere + [2, 0, 0], "Wavenumber", 2.0);
        details.responseNorm = norm(op.apply(ones(size(op.matrix, 2), 1)));
        failures = requirePass(failures, isfinite(details.responseNorm) && details.responseNorm > 0, ...
            "Mid-ka pulsating sphere toy response is invalid.");
    case 5
        K = HelmholtzKernel(target, source, "Wavenumber", 0.05, ...
            "SourceNormals", repmat([1, 0, 0], size(source, 1), 1));
        details.doubleLayerNorm = norm(K.doubleLayerSourceNormal, "fro");
        failures = requirePass(failures, isfinite(details.doubleLayerNorm), "Rigid-sphere small-ka toy derivative is invalid.");
    case 6
        op = AcousticSingleLayer(target, source, "Wavenumber", 1.5);
        rhs = exp(1i * target(:, 1));
        density = op.matrix \ rhs;
        details.solutionNorm = norm(density);
        failures = requirePass(failures, all(isfinite(density)), "Dirichlet dense acoustic solve produced invalid density.");
    case 7
        ks = [0.5, 1.0, 1.5];
        phases = zeros(size(ks));
        for k = 1:numel(ks)
            K = HelmholtzKernel([0, 0, 0], [1, 0, 0], "Wavenumber", ks(k));
            phases(k) = angle(K.singleLayer(1, 1));
        end
        details.phases = phases;
        failures = requirePass(failures, all(diff(phases) > 0), "Frequency sweep phase did not increase.");
    case 8
        op = AcousticSingleLayer(target, source, "Wavenumber", 2.5);
        x = sin((1:size(source, 1)).');
        details.reproducibilityError = norm(op.apply(x) - op.matrix * x);
        failures = requirePass(failures, details.reproducibilityError < tol, "Complex acoustic matvec is not reproducible.");
    case 9
        K = HelmholtzKernel(target, source, "Wavenumber", 1.0);
        details.matrixNorm = norm(K.singleLayer, "fro");
        failures = requirePass(failures, isfinite(details.matrixNorm) && details.matrixNorm > 0, ...
            "Helmholtz H-matrix candidate matrix is invalid.");
    case 10
        op = AcousticSingleLayer(target, source, "Wavenumber", 1.0);
        details.referenceNorm = norm(op.matrix, "fro");
        failures = requirePass(failures, caps.hasBem && caps.hasHelmholtzSL && details.referenceNorm > 0, ...
            "NGSolve Helmholtz BEM comparison gate is unavailable.");
end
end


function [failures, details] = verifyScalarCoupling(item, caps, failures)
model = unitModel();
model = model.assemble();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "scalar_fem_bem_coupling", "slot", slot);
switch slot
    case 1
        T = model.operators.trace.matrix;
        details.traceSize = size(T);
        failures = requirePass(failures, isequal(size(T), [4, 4]), "Unit tetra trace matrix size is wrong.");
    case 2
        patch = fourTetModel();
        patch = patch.assemble();
        traceNodes = patch.operators.trace.femNodeIds;
        allNodes = (1:size(patch.mesh.vtx, 1)).';
        interior = setdiff(allNodes, traceNodes);
        details.interiorNodes = interior.';
        failures = requirePass(failures, ~isempty(interior), "Interior node was not separated from the trace.");
    case 3
        u = sum(model.mesh.vtx, 2);
        residual = model.operators.fem.stiffness * u;
        traceResidual = model.operators.trace.matrix * residual;
        details.traceResidualNorm = norm(traceResidual);
        failures = requirePass(failures, all(isfinite(traceResidual)), "Dirichlet-to-Neumann toy residual is invalid.");
    case 4
        surfaceMass = model.operators.bem.surfaceMass;
        details.surfaceMassTrace = trace(surfaceMass);
        failures = requirePass(failures, issymmetricWithin(surfaceMass, tol), "Boundary P1 surface mass is not symmetric.");
    case 5
        balance = norm(full(model.operators.fem.stiffness) * ones(size(model.mesh.vtx, 1), 1));
        details.balance = balance;
        failures = requirePass(failures, balance < tol, "FEM/BEM flux balance for constants failed.");
    case 6
        [~, scaled] = model.h1.stiffness(2);
        err = norm(full(scaled.stiffness - 2 * model.operators.fem.stiffness), "fro");
        details.error = err;
        failures = requirePass(failures, err < tol, "Single-material coupling stiffness scaling failed.");
    case 7
        patch = fourTetModel();
        coeff = 1:size(patch.mesh.tet, 1);
        [~, K] = patch.h1.stiffness(coeff);
        details.nonzeros = nnz(K.stiffness);
        failures = requirePass(failures, issymmetricWithin(K.stiffness, tol), "Two-material coupling stiffness is not symmetric.");
    case 8
        sphere = spherePoints(12);
        K = HelmholtzKernel(sphere, sphere + [2, 0, 0]);
        details.meanPotential = mean(K.singleLayer * ones(size(sphere, 1), 1));
        failures = requirePass(failures, details.meanPotential > 0, "Exterior scalar potential toy is invalid.");
    case 9
        [target, source] = separatedPointClouds();
        H = HMatrix(target, source);
        y = H.matvec(ones(size(source, 1), 1));
        details.responseNorm = norm(y);
        failures = requirePass(failures, details.responseNorm > 0, "Scalar FEM/BEM H-matrix response is invalid.");
    case 10
        details.h1Dofs = caps.h1Dofs;
        details.hasLaplaceSL = caps.hasLaplaceSL;
        failures = requirePass(failures, caps.h1Dofs == 4 && caps.hasLaplaceSL, ...
            "NGSolve coupled scalar comparison gate is unavailable.");
end
end


function [failures, details] = verifyRwgHCurlTrace(item, caps, failures)
model = unitModel();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "rwg_hcurl_trace", "slot", slot);
switch slot
    case 1
        details.rwgDofs = numel(model.rwg.dofEdgeIds);
        failures = requirePass(failures, numel(model.rwg.dofEdgeIds) == 6, ...
            "Closed tetrahedron should have six RWG dofs.");
    case 2
        openSurface = model.surface;
        openSurface.tri = openSurface.tri(1, :);
        openRwg = RwgSpace(openSurface);
        details.openRwgDofs = numel(openRwg.dofEdgeIds);
        failures = requirePass(failures, isempty(openRwg.dofEdgeIds), ...
            "Open triangle should not create interior RWG edge dofs.");
    case 3
        details.oppositeVertices = model.rwg.oppositeVerticesLocal;
        failures = requirePass(failures, all(model.rwg.oppositeVerticesLocal(:) > 0), ...
            "RWG opposite-vertex map is incomplete.");
    case 4
        ids = model.rwgToHcurlEdgeIds;
        details.hcurlEdgeIds = ids.';
        failures = requirePass(failures, all(ids >= 1 & ids <= size(model.hcurl.edges, 1)), ...
            "RWG to HCurl edge ids are out of range.");
    case 5
        signs = model.rwg.triEdgeSigns;
        details.signs = signs;
        failures = requirePass(failures, all(abs(signs(:)) == 1), "Boundary edge orientation signs must be +/-1.");
    case 6
        patch = fourTetModel();
        details.boundaryDofs = numel(patch.rwg.dofEdgeIds);
        failures = requirePass(failures, numel(patch.rwg.dofEdgeIds) > 0, ...
            "Two-tet boundary should expose RWG dofs.");
    case 7
        current = ones(numel(model.rwg.dofEdgeIds), 1);
        lifted = zeros(size(model.hcurl.edges, 1), 1);
        lifted(model.rwgToHcurlEdgeIds) = current;
        details.currentNorm = norm(lifted);
        failures = requirePass(failures, norm(lifted) > 0, "Surface current toy did not lift to HCurl edges.");
    case 8
        [~, ~, ned] = model.hcurl.matrices();
        traceIds = model.rwgToHcurlEdgeIds;
        traceMass = ned.mass(traceIds, traceIds);
        details.traceMassNorm = norm(full(traceMass), "fro");
        failures = requirePass(failures, details.traceMassNorm > 0, "Magnetic vector trace toy is invalid.");
    case 9
        edgeTrace = sparse(1:numel(model.rwgToHcurlEdgeIds), model.rwgToHcurlEdgeIds, ...
            1, numel(model.rwgToHcurlEdgeIds), size(model.hcurl.edges, 1));
        details.edgeTraceSize = size(edgeTrace);
        failures = requirePass(failures, nnz(edgeTrace) == numel(model.rwgToHcurlEdgeIds), ...
            "Edge trace sparse map has wrong nonzero count.");
    case 10
        [~, ~, ned] = model.hcurl.matrices();
        err = norm(full(ned.mass) - caps.hcurlMass, "fro");
        details.hcurlMassError = err;
        failures = requirePass(failures, caps.hcurlDofs == 6 && err < tol, ...
            "NGSolve HCurl trace comparison gate failed.");
end
end


function [failures, details] = verifyNgsolveBemReference(item, caps, failures)
model = unitModel();
slot = caseSlot(item);
tol = item.tolerance;

details = struct("gate", "ngsolve_bem_reference", "slot", slot);
switch slot
    case 1
        details.meshVertices = caps.meshVertices;
        details.meshElements = caps.meshElements;
        failures = requirePass(failures, caps.meshVertices == 4 && caps.meshElements == 1 && caps.meshEdges == 6, ...
            "NGSolve mesh import smoke failed.");
    case 2
        [~, K] = model.h1.stiffness();
        err = norm(full(K.stiffness) - caps.h1Stiffness, "fro");
        details.error = err;
        failures = requirePass(failures, caps.h1Dofs == 4 && err < tol, "NGSolve P1 scalar smoke failed.");
    case 3
        [~, ~, N] = model.hcurl.matrices();
        err = norm(full(N.curlCurl) - caps.hcurlCurlCurl, "fro");
        details.error = err;
        failures = requirePass(failures, caps.hcurlDofs == 6 && err < tol, "NGSolve HCurl smoke failed.");
    case 4
        details.hasLaplaceSL = caps.hasLaplaceSL;
        failures = requirePass(failures, caps.hasBem && caps.hasLaplaceSL, "NGSolve Laplace BEM smoke failed.");
    case 5
        details.hasHelmholtzSL = caps.hasHelmholtzSL;
        failures = requirePass(failures, caps.hasBem && caps.hasHelmholtzSL, "NGSolve Helmholtz BEM smoke failed.");
    case 6
        K = HelmholtzKernel([0, 0, 0], [1, 0, 0], "Wavenumber", 1e-8);
        details.lowFrequencyValue = K.singleLayer(1, 1);
        failures = requirePass(failures, caps.hasHelmholtzSL && isfinite(K.singleLayer(1, 1)), ...
            "NGSolve low-frequency smoke gate failed.");
    case 7
        sphere = spherePoints(16);
        K = HelmholtzKernel(sphere, sphere + [2, 0, 0]);
        details.meanPotential = mean(K.singleLayer * ones(size(sphere, 1), 1));
        failures = requirePass(failures, details.meanPotential > 0 && caps.hasLaplaceSL, ...
            "Sphere analytic Laplace gate failed.");
    case 8
        sphere = spherePoints(16);
        K = HelmholtzKernel(sphere, sphere + [2, 0, 0], "Wavenumber", 1.0);
        details.meanAmplitude = mean(abs(K.singleLayer * ones(size(sphere, 1), 1)));
        failures = requirePass(failures, details.meanAmplitude > 0 && caps.hasHelmholtzSL, ...
            "Sphere analytic acoustic gate failed.");
    case 9
        coupled = model.assemble();
        details.traceRows = size(coupled.operators.trace.matrix, 1);
        failures = requirePass(failures, size(coupled.operators.trace.matrix, 1) == 4 && caps.hasLaplaceSL, ...
            "NGSolve FEM/BEM coupling smoke gate failed.");
    case 10
        [~, K] = model.h1.stiffness();
        [~, M] = model.h1.mass();
        [~, ~, N] = model.hcurl.matrices();
        details.h1Error = norm(full(K.stiffness) - caps.h1Stiffness, "fro") + ...
            norm(full(M.mass) - caps.h1Mass, "fro");
        details.hcurlError = norm(full(N.mass) - caps.hcurlMass, "fro") + ...
            norm(full(N.curlCurl) - caps.hcurlCurlCurl, "fro");
        failures = requirePass(failures, details.h1Error < tol && details.hcurlError < tol && ...
            caps.hasLaplaceSL && caps.hasHelmholtzSL, "Full pipeline capstone failed.");
end
end


function failures = requirePass(failures, condition, message)
if ~condition
    failures(end + 1, 1) = string(message);
end
end


function slot = caseSlot(item)
number = sscanf(char(item.id), "GYP-%d");
slot = mod(number - 1, 10) + 1;
end


function tf = issymmetricWithin(A, tol)
tf = norm(full(A - A.'), "fro") < tol;
end


function model = unitModel()
model = FemBemModel(fullfile(gypsilabRepoRoot(), "fixtures", "mesh_topology", "unit_tetra.vol"));
end


function model = fourTetModel()
model = FemBemModel(fullfile(gypsilabRepoRoot(), "fixtures", "mesh_topology", "four_tet_interior_node.vol"));
end


function area = triangleArea(points)
area = 0.5 * norm(cross(points(2, :) - points(1, :), points(3, :) - points(1, :)));
end


function [centers, areas] = triangleCentroidsAndAreas(vtx, tri)
a = vtx(tri(:, 1), :);
b = vtx(tri(:, 2), :);
c = vtx(tri(:, 3), :);
centers = (a + b + c) / 3;
areas = 0.5 * sqrt(sum(cross(b - a, c - a, 2).^2, 2));
end


function G = edgeGradientMatrix(edges, nNodes)
nEdges = size(edges, 1);
ii = repelem((1:nEdges).', 2);
jj = reshape(edges.', [], 1);
vv = repmat([-1; 1], nEdges, 1);
G = sparse(ii, jj, vv, nEdges, nNodes);
end


function [target, source] = separatedPointClouds()
[x, y] = ndgrid(linspace(0, 1, 4), linspace(0, 1, 4));
source = [x(:), y(:), zeros(numel(x), 1)];
target = source + [3, 0.25, 0.5];
end


function d = pairwiseDistance(target, source)
d = zeros(size(target, 1), size(source, 1));
for i = 1:size(target, 1)
    delta = source - target(i, :);
    d(i, :) = sqrt(sum(delta.^2, 2)).';
end
end


function points = spherePoints(n)
theta = linspace(0, 2 * pi, n + 1).';
theta(end) = [];
z = linspace(-0.75, 0.75, n).';
r = sqrt(max(0, 1 - z.^2));
points = [r .* cos(theta), r .* sin(theta), z];
end


function points = cubeSurfacePoints()
base = [-1, -1; 1, -1; -1, 1; 1, 1];
points = [
    base, -ones(4, 1)
    base, ones(4, 1)
    base(:, 1), -ones(4, 1), base(:, 2)
    base(:, 1), ones(4, 1), base(:, 2)
    -ones(4, 1), base
    ones(4, 1), base
];
end
