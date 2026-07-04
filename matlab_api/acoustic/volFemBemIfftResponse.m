function result = volFemBemIfftResponse(volFile, options)
%VOLFEMBEMIFFTRESPONSE Time response from .vol P1 FEM/BEM frequency solves.
%
%   result = volFemBemIfftResponse("model.vol") builds a FemBemModel from a
%   first-order Netgen .vol mesh, solves the scalar acoustic FEM/BEM
%   transmission problem on a frequency grid, multiplies the transfer
%   function by a real Ricker-pulse spectrum, and reconstructs exterior
%   observation pressure with an inverse FFT.
%
%   This is the readable "real solver" bridge after reduced drum cartoons:
%   .vol volume tetrahedra -> H1/P1 FEM, .vol boundary triangles -> P1 BEM,
%   frequency-domain Helmholtz FEM/BEM -> time-domain response by iFFT.  It
%   is not convolution-quadrature TD-BEM yet; it is the frequency-sweep/iFFT
%   route, deliberately kept short enough for students to read.

arguments
    volFile (1,1) string = ""
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 7)} = 32
    options.FinalTime (1,1) double {mustBePositive} = 20.0
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.InteriorSoundSpeed double = []
    options.DensityRatio (1,1) double {mustBePositive} = 1.0
    options.VolumeSource (1,1) double = 1.0
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
    options.ObservationPoints double = []
    options.PulseCenterTime double = []
    options.PulseWidth double = []
end

if mod(options.NumTime, 2) ~= 0
    error("volFemBemIfftResponse:NumTime", ...
        "NumTime must be even so the positive and negative frequency bins pair cleanly.");
end

if strlength(volFile) == 0
    volFile = defaultFixture();
end

meshTimer = tic;
model = FemBemModel(volFile);
meshSeconds = toc(meshTimer);

if isempty(options.ObservationPoints)
    obs = defaultObservationPoints(model.mesh.vtx);
else
    obs = options.ObservationPoints;
end
if size(obs, 2) ~= 3
    error("volFemBemIfftResponse:ObservationPoints", ...
        "ObservationPoints must have three columns.");
end

N = options.NumTime;
T = options.FinalTime;
dt = T / N;
t = (0:N-1).' * dt;
source = rickerPulse(t, T, options.PulseCenterTime, options.PulseWidth);
sourceSpectrum = fft(source);

nHalf = N / 2 + 1;
positiveBins = 1:(nHalf - 1);          % DC through last non-Nyquist bin
omega = 2 * pi * (positiveBins - 1).' / T;
k0 = omega / options.SoundSpeed;
if isempty(options.InteriorSoundSpeed)
    k1 = k0;
else
    k1 = omega / options.InteriorSoundSpeed;
end

solveTimer = tic;
H = zeros(numel(positiveBins), size(obs, 1));
residuals = zeros(numel(positiveBins), 1);
status = strings(numel(positiveBins), 1);
for b = 1:numel(positiveBins)
    if k0(b) == 0
        sol = femBemCoupledSolve(model, ...
            "VolumeSource", options.VolumeSource, ...
            "QuadratureOrder", options.QuadratureOrder);
    else
        sol = femBemCoupledSolve(model, ...
            "VolumeSource", options.VolumeSource, ...
            "Wavenumber", k0(b), ...
            "InteriorWavenumber", k1(b), ...
            "DensityRatio", options.DensityRatio, ...
            "QuadratureOrder", options.QuadratureOrder);
    end
    H(b, :) = sol.exteriorPotentialAt(obs).';
    residuals(b) = sol.residualNorm;
    status(b) = sol.status;
end
solveSeconds = toc(solveTimer);

postTimer = tic;
pressureSpectrum = zeros(N, size(obs, 1));
pressureSpectrum(positiveBins, :) = H .* sourceSpectrum(positiveBins);
for b = 2:numel(positiveBins)
    pressureSpectrum(N - b + 2, :) = conj(pressureSpectrum(b, :));
end
% The Nyquist bin is self-conjugate for a real time series.  Keeping it zero
% avoids pretending that the one-sided Helmholtz solve supplies both signs.
nyquistBin = nHalf;
pressureSpectrum(nyquistBin, :) = 0;
timeComplex = ifft(pressureSpectrum, [], 1);
pressure = real(timeComplex);
postSeconds = toc(postTimer);

result = struct();
result.kind = "vol_p1_fem_bem_frequency_ifft_time_response";
result.method = "frequency_domain_p1_fem_p1_bem_plus_inverse_fft";
result.volFile = string(volFile);
result.meshId = model.mesh.meshId;
result.meshSummary = model.mesh.summary;
result.meshSourceId = model.mesh.sourceFileId;
result.policy = "first_order_vol_tri_tet_p1_volume_fem_p1_surface_bem";
result.time = t;
result.timeStep = dt;
result.finalTime = T;
result.sourceTime = source;
result.sourceSpectrumPositive = sourceSpectrum(positiveBins);
result.observationPoints = obs;
result.frequency = omega / (2 * pi);
result.angularFrequency = omega;
result.wavenumber = k0;
result.interiorWavenumber = k1;
result.frequencyResponse = H;
result.pressureSpectrum = pressureSpectrum;
result.pressure = pressure;
result.residuals = residuals;
result.frequencySolveStatus = status;
result.timing = struct( ...
    "mesh_or_import", meshSeconds, ...
    "frequency_solves", solveSeconds, ...
    "ifft_postprocess", postSeconds);
result.summary = struct( ...
    "num_time", N, ...
    "num_frequency_solves", numel(positiveBins), ...
    "num_observation_points", size(obs, 1), ...
    "max_abs_pressure", max(abs(pressure), [], "all"), ...
    "max_solve_residual", max(residuals), ...
    "max_imag_time_pressure_before_real", max(abs(imag(timeComplex)), [], "all"));
result.checks = struct( ...
    "vol_mesh_tri_tet", model.mesh.summary.triangles > 0 && model.mesh.summary.tets > 0, ...
    "p1_volume_fem", model.h1.order == 1 && model.h1.cell == "tetrahedron", ...
    "p1_boundary_bem", model.scalarBem.order == 1 && model.scalarBem.cell == "triangle", ...
    "frequency_solves_ok", all(status == "ok"), ...
    "finite_pressure", all(isfinite(pressure), "all"), ...
    "nonzero_time_response", max(abs(pressure), [], "all") > 0, ...
    "ifft_is_real", result.summary.max_imag_time_pressure_before_real < 1e-10 * max(1, result.summary.max_abs_pressure), ...
    "nyquist_bin_zeroed", all(pressureSpectrum(nyquistBin, :) == 0));
if all(structfun(@(v) logical(v), result.checks))
    result.status = "ok";
else
    result.status = "needs_attention";
end
end


function volFile = defaultFixture()
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_tetra.vol"));
end


function obs = defaultObservationPoints(nodes)
center = mean(nodes, 1);
r = max(vecnorm(nodes - center, 2, 2));
if r <= 0
    r = 1;
end
obs = center + r * [
    0.0 0.0 2.5
    2.5 0.0 0.0
    0.0 2.5 1.5
    ];
end


function source = rickerPulse(t, finalTime, centerTime, width)
if isempty(centerTime)
    centerTime = 0.22 * finalTime;
end
if isempty(width)
    width = 0.055 * finalTime;
end
x = (t - centerTime) / width;
source = (1 - 2 * x.^2) .* exp(-x.^2);
source = source - mean(source);
end
