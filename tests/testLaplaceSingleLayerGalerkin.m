function tests = testLaplaceSingleLayerGalerkin
%TESTLAPLACESINGLELAYERGALERKIN Laplace single layer vs Helmholtz on imag axis.
%
% The Lubich CQ time-domain BEM (volTdBemConvolutionQuadrature) is only as
% trustworthy as its Laplace-domain building block V(s).  This locks that block
% to the analytically-validated frequency-domain single layer: on the imaginary
% axis s = -1i c k the retarded kernel exp(-s r/c) equals the Helmholtz kernel
% exp(1i k r), so V(-1i c k) equals GalerkinSingleLayer(k) plus the KNOWN
% coincident-quadrature-node limit term Delta (this operator keeps the finite
% smooth-correction limit; GalerkinSingleLayer drops it).  Matching to machine
% precision pins the s/c scaling, the exponent sign, and the 1/(4 pi) - not just
% a self-consistent CQ residual.

tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function surface = coarseSphere()
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_coarse.vol");
surface = VolMesh(volFile).boundary();
end


function testImaginaryAxisMatchesHelmholtzPlusCoincidentLimit(testCase)
% s = -1i c k  =>  exp(-s r/c) = exp(1i k r).  V(s) must equal the Helmholtz
% single layer plus the coincident-node limit term Delta_ij =
% (-alpha/4pi) sum_g w_g^2 phi_i(x_g) phi_j(x_g) with alpha = s/c (= -1i k here,
% so -alpha = 1i k), to machine precision, at every quadrature order and k.
surface = coarseSphere();
c = 1.3;
for q = [1 3 7]
    for k = [0.5 1.5 3.0]
        s = -1i * c * k;
        V = laplaceSingleLayerGalerkin(surface, s, c, q);
        Vref = GalerkinSingleLayer(surface, "Wavenumber", k, ...
            "QuadratureOrder", q).matrix;
        quad = SurfaceQuadrature(surface, q);
        w = quad.weights;
        B = quad.basis;
        Delta = (1i * k) / (4 * pi) * (B.' * (w.^2 .* B));
        relErr = norm(V - (Vref + Delta), "fro") / norm(Vref, "fro");
        verifyLessThan(testCase, relErr, 1e-12, ...
            sprintf("V(-i c k) != Helmholtz + Delta at q=%d k=%.2f (relerr %.2e)", ...
            q, k, relErr));
    end
end
end


function testOperatorIsGenuinelyFrequencyDependent(testCase)
% Guard against a trivial pass where the smooth correction collapses and V(s)
% is just the real k=0 Laplace single layer: the correction must be a nonzero
% complex contribution that grows with |s|.
surface = coarseSphere();
c = 1.0;
V0 = laplaceSingleLayerGalerkin(surface, 0, c, 3);            % Laplace limit
verifyLessThan(testCase, norm(imag(V0), "fro"), 1e-12);       % real at s = 0
verifyGreaterThan(testCase, norm(V0, "fro"), 0);
Vlow  = laplaceSingleLayerGalerkin(surface, -1i * 0.5, c, 3);
Vhigh = laplaceSingleLayerGalerkin(surface, -1i * 3.0, c, 3);
corrLow  = norm(Vlow  - V0, "fro") / norm(V0, "fro");
corrHigh = norm(Vhigh - V0, "fro") / norm(V0, "fro");
verifyGreaterThan(testCase, corrLow, 1e-3);                   % genuinely complex
verifyGreaterThan(testCase, corrHigh, corrLow);              % grows with |s|
end


function testSoundSpeedEntersOnlyAsRatio(testCase)
% The retarded kernel exp(-s r/c) depends on s and c only through s/c, so
% V(s; c) = V(s/c; 1) for a general complex s (a missing /c would break this).
surface = coarseSphere();
s = 0.7 - 2.1i;                       % positive real part, like a real CQ node
Va = laplaceSingleLayerGalerkin(surface, s, 2.0, 3);
Vb = laplaceSingleLayerGalerkin(surface, s / 2.0, 1.0, 3);
verifyLessThan(testCase, norm(Va - Vb, "fro") / norm(Vb, "fro"), 1e-12);
end


function testNonFiniteLaplaceParameterFailsLoud(testCase)
% The promoted public operator rejects a non-finite s rather than silently
% returning NaN/Inf (fail-loud discipline).
surface = coarseSphere();
verifyError(testCase, @() laplaceSingleLayerGalerkin(surface, complex(Inf, 0), 1.0, 1), ...
    "MATLAB:validators:mustBeFinite");
verifyError(testCase, @() laplaceSingleLayerGalerkin(surface, complex(0, NaN), 1.0, 1), ...
    "MATLAB:validators:mustBeFinite");
end


function testShippedCqSolveUsesThisOperator(testCase)
% End-to-end: rebuild the CQ boundary density from laplaceSingleLayerGalerkin at
% the CQ's OWN Laplace nodes and confirm it reproduces the shipped
% result.boundaryDensity.  This protects the operator the CQ actually uses - a
% mutated operator (or a reverted private copy in volTdBem) would diverge here.
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol");
N = 8;
q = 3;
result = volTdBemConvolutionQuadrature(volFile, "NumTime", N, "TimeStep", 0.6, ...
    "QuadratureOrder", q, "Method", "BDF1");
surface = VolMesh(volFile).boundary();
nB = size(surface.vtx, 1);
n = (0:N-1).';
rho = result.cqRadius;
s = result.cqLaplaceParameter;
boundaryHat = fft((rho .^ n) .* result.boundaryData, [], 1);
densityHat = zeros(N, nB);
for ell = 1:N
    V = laplaceSingleLayerGalerkin(surface, s(ell), 1.0, q);   % c = 1 (CQ default)
    densityHat(ell, :) = (V \ boundaryHat(ell, :).').';
end
density = real((rho .^ (-n)) .* ifft(densityHat, [], 1));
relErr = norm(density - result.boundaryDensity, "fro") / ...
    max(1, norm(result.boundaryDensity, "fro"));
verifyLessThan(testCase, relErr, 1e-10);
verifyEqual(testCase, result.status, "ok");
end
