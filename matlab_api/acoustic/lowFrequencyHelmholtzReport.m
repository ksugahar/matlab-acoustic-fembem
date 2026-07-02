function report = lowFrequencyHelmholtzReport(targetPoint, sourcePoint, options)
%LOWFREQUENCYHELMHOLTZREPORT Explain the low-frequency BEM split.
%
% This tiny report is intentionally more verbose than a production kernel.
% Students can see the singular Laplace term, the smooth Helmholtz correction,
% and the cancellation scale that makes direct subtraction fragile when k*r is
% small.

arguments
    targetPoint (1,3) double
    sourcePoint (1,3) double
    options.Wavenumber (1,1) double {mustBeNonnegative} = 1e-6
    options.TaylorCutoff (1,1) double {mustBePositive} = 1e-4
    options.SeriesTerms (1,1) double {mustBeInteger, mustBePositive} = 8
end

delta = targetPoint - sourcePoint;
r = norm(delta);
if r <= 0
    error("lowFrequencyHelmholtzReport:distance", ...
        "targetPoint and sourcePoint must be distinct.");
end

k = options.Wavenumber;
kernel = HelmholtzKernel(targetPoint, sourcePoint, ...
    "Wavenumber", k, ...
    "TaylorCutoff", options.TaylorCutoff, ...
    "SeriesTerms", options.SeriesTerms);

directGreen = exp(1i * k * r) / (4 * pi * r);
laplaceTerm = kernel.singleLayerLaplace;
stableCorrection = kernel.singleLayerCorrection;
directCorrection = directGreen - laplaceTerm;

report = struct();
report.kind = "low_frequency_helmholtz_teaching_report";
report.policy = "readable_bem_kernel_split_not_production_quadrature";
report.distance = r;
report.wavenumber = k;
report.kr = k * r;
report.laplaceTerm = laplaceTerm;
report.stableCorrection = stableCorrection;
report.directCorrection = directCorrection;
report.singleLayer = kernel.singleLayer;
report.directGreen = directGreen;
report.stableError = abs(kernel.singleLayer - directGreen);
report.correctionAgreement = abs(stableCorrection - directCorrection);
report.cancellationRatio = abs(laplaceTerm) / max(abs(stableCorrection), eps);
report.notes = [
    "G_k = G_0 + (exp(1i*k*r)-1)/(4*pi*r)"
    "G_0 keeps the singular Laplace quadrature visible"
    "the correction is smooth and evaluated by Taylor/expm1 for low k*r"
    ];
end
