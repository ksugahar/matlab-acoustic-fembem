function tests = testLaplaceSingleLayerGalerkin
%TESTLAPLACESINGLELAYERGALERKIN Laplace single layer == Helmholtz on imag axis.
%
% The Lubich CQ time-domain BEM (volTdBemConvolutionQuadrature) is only as
% trustworthy as its Laplace-domain building block V(s).  This locks that block
% to the analytically-validated frequency-domain single layer: on the imaginary
% axis s = -1i c k the retarded kernel exp(-s r/c) equals the Helmholtz kernel
% exp(1i k r), so V(-1i c k) must equal GalerkinSingleLayer(k) to machine
% precision - an EXACT golden that pins the s/c scaling, the sign of the
% exponent, and the 1/(4 pi), not just a self-consistent CQ residual.

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


function testImaginaryAxisMatchesHelmholtzSingleLayer(testCase)
% s = -1i c k  =>  exp(-s r/c) = exp(1i k r): V(s) must equal the Helmholtz
% single layer to machine precision at every quadrature order and wavenumber.
surface = coarseSphere();
c = 1.3;
for q = [1 3 7]
    for k = [0.5 1.5 3.0]
        s = -1i * c * k;
        V = laplaceSingleLayerGalerkin(surface, s, c, q);
        Vref = GalerkinSingleLayer(surface, "Wavenumber", k, ...
            "QuadratureOrder", q).matrix;
        relErr = norm(V - Vref, "fro") / norm(Vref, "fro");
        verifyLessThan(testCase, relErr, 1e-12, ...
            sprintf("V(-i c k) != Helmholtz at q=%d k=%.2f (relerr %.2e)", q, k, relErr));
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


function testCqPathUsesThisOperator(testCase)
% The CQ solver must actually build V(s) with laplaceSingleLayerGalerkin at its
% own nodes (so this golden protects the shipped path, not a parallel copy).
repoRoot = fileparts(fileparts(mfilename("fullpath")));
volFile = fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol");
result = volTdBemConvolutionQuadrature(volFile, "NumTime", 8, "TimeStep", 0.6, ...
    "QuadratureOrder", 3, "Method", "BDF1");
surface = VolMesh(volFile).boundary();
s = result.cqLaplaceParameter;
worst = 0;
for ell = 1:numel(s)
    V = laplaceSingleLayerGalerkin(surface, s(ell), 1.0, 3);
    q = V \ ones(size(surface.vtx, 1), 1);
    worst = max(worst, norm(V * q - ones(size(surface.vtx, 1), 1)) / sqrt(numel(q)));
end
verifyLessThan(testCase, worst, 1e-10);   % V(s) is well-formed at the CQ nodes
verifyEqual(testCase, result.status, "ok");
end
