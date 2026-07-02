function report = educationalLowFrequencyHelmholtzKernelManifest(targetPoint, sourcePoint, options)
%EDUCATIONALLOWFREQUENCYHELMHOLTZKERNELMANIFEST Manifest for low-frequency BEM kernels.
%
% This wrapper keeps the readable identity of the split visible:
% Helmholtz single-layer = Laplace singular kernel + smooth Taylor/expm1
% correction.  It is meant for teaching notebooks and validation manifests,
% not production quadrature.

arguments
    targetPoint (1,3) double
    sourcePoint (1,3) double
    options.Wavenumber (1,1) double {mustBeNonnegative} = 1e-9
    options.KernelFamily (1,1) string = "helmholtz_single_layer"
    options.ExpectedKernelFamily (1,1) string = "helmholtz_single_layer"
    options.LowFrequencyStrategy (1,1) string = "laplace_plus_expm1_taylor_correction"
    options.ExpectedLowFrequencyStrategy (1,1) string = "laplace_plus_expm1_taylor_correction"
    options.TimeConvention (1,1) string = "exp(+i*k*r) MATLAB teaching convention"
    options.ExpectedTimeConvention (1,1) string = "exp(+i*k*r) MATLAB teaching convention"
    options.MaxKr (1,1) double {mustBePositive} = 1e-6
    options.MaxStableError (1,1) double {mustBeNonnegative} = 1e-12
    options.MaxCorrectionAgreement (1,1) double {mustBeNonnegative} = 1e-12
end

base = lowFrequencyHelmholtzTeachingReport(targetPoint, sourcePoint, ...
    "Wavenumber", options.Wavenumber);

report = base;
report.kind = "low_frequency_helmholtz_kernel_manifest";
report.policy = "readable_low_frequency_helmholtz_kernel_manifest_gate";
report.kernelFamily = options.KernelFamily;
report.expectedKernelFamily = options.ExpectedKernelFamily;
report.lowFrequencyStrategy = options.LowFrequencyStrategy;
report.expectedLowFrequencyStrategy = options.ExpectedLowFrequencyStrategy;
report.timeConvention = options.TimeConvention;
report.expectedTimeConvention = options.ExpectedTimeConvention;
report.krAbs = abs(base.kr);
report.maxKr = options.MaxKr;
report.maxStableError = options.MaxStableError;
report.maxCorrectionAgreement = options.MaxCorrectionAgreement;

checks = struct();
checks.kernelFamilyRecorded = strlength(report.kernelFamily) > 0;
checks.kernelFamilyMatchesExpected = report.kernelFamily == report.expectedKernelFamily;
checks.lowFrequencyStrategyRecorded = strlength(report.lowFrequencyStrategy) > 0;
checks.lowFrequencyStrategyMatchesExpected = report.lowFrequencyStrategy == report.expectedLowFrequencyStrategy;
checks.timeConventionRecorded = strlength(report.timeConvention) > 0;
checks.timeConventionMatchesExpected = report.timeConvention == report.expectedTimeConvention;
checks.krWithinLowFrequencyLimit = report.krAbs <= report.maxKr;
checks.laplaceTermRecorded = ~isempty(report.laplaceTerm);
checks.stableCorrectionRecorded = ~isempty(report.stableCorrection);
checks.directCorrectionRecorded = ~isempty(report.directCorrection);
checks.stableErrorWithinTolerance = report.stableError <= report.maxStableError;
checks.correctionAgreementWithinTolerance = report.correctionAgreement <= report.maxCorrectionAgreement;
report.checks = checks;

flags = structfun(@(value) isequal(value, true), checks);
if all(flags)
    report.status = "ok";
else
    report.status = "needs_attention";
end
report.notes = [
    report.notes
    "kernelFamily, timeConvention, and lowFrequencyStrategy are part of the value identity"
    "krAbs must remain inside the low-frequency teaching range before this split is used as evidence"
    ];
end
