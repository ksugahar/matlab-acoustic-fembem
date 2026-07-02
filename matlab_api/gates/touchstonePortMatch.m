function result = touchstonePortMatch(s11, options)
%touchstonePortMatch Readable S11 match-quality teaching gate.
%
% This helper keeps one-port match metrics separate from two-port insertion
% loss.  It is useful before CST/Touchstone, ngsolve.bem, or measurement rows
% are ranked by match quality or reused in optimization notebooks.
%
%   |Gamma| = |S11|
%   VSWR = (1 + |Gamma|) / (1 - |Gamma|)
%   return_loss_db = -20*log10(|Gamma|)
%   mismatch_loss_db = -10*log10(1 - |Gamma|^2)

arguments
    s11 (1,1) double
    options.ReturnLossMinDb (1,1) double = NaN
    options.VswrMax (1,1) double = NaN
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-12
end

gamma = abs(s11);
tolerance = options.Tolerance;
passiveReflectionOk = gamma <= 1 + tolerance;

if gamma >= 1
    vswr = Inf;
else
    vswr = (1 + gamma) / (1 - gamma);
end

if gamma == 0
    returnLossDb = Inf;
else
    returnLossDb = -20 * log10(gamma);
end

transmittedPowerFraction = 1 - gamma^2;
reflectedPowerFraction = gamma^2;
if transmittedPowerFraction > 0
    mismatchLossDb = -10 * log10(transmittedPowerFraction);
else
    mismatchLossDb = Inf;
end

if isinf(vswr)
    gammaFromVswr = 1;
else
    gammaFromVswr = (vswr - 1) / (vswr + 1);
end
if isinf(returnLossDb)
    gammaFromReturnLoss = 0;
else
    gammaFromReturnLoss = 10^(-returnLossDb / 20);
end

returnLossLimitOk = true;
if ~isnan(options.ReturnLossMinDb)
    returnLossLimitOk = returnLossDb >= options.ReturnLossMinDb;
end
vswrLimitOk = true;
if ~isnan(options.VswrMax)
    vswrLimitOk = vswr <= options.VswrMax;
end

checks = struct( ...
    "passiveReflectionOk", passiveReflectionOk, ...
    "vswrRoundTripOk", abs(gammaFromVswr - min(gamma, 1)) <= tolerance, ...
    "returnLossRoundTripOk", abs(gammaFromReturnLoss - gamma) <= tolerance, ...
    "transmittedPowerFractionNonnegative", transmittedPowerFraction >= -tolerance, ...
    "returnLossLimitOk", returnLossLimitOk, ...
    "vswrLimitOk", vswrLimitOk);

result = struct();
result.kind = "touchstone_port_match";
result.policy = "readable_one_port_match_quality_gate";
result.s11 = s11;
result.reflectionCoefficient = gamma;
result.vswr = vswr;
result.returnLossDb = returnLossDb;
result.mismatchLossDb = mismatchLossDb;
result.transmittedPowerFraction = transmittedPowerFraction;
result.reflectedPowerFraction = reflectedPowerFraction;
result.returnLossMinDb = options.ReturnLossMinDb;
result.vswrMax = options.VswrMax;
result.checks = checks;
result.notes = [
    "S11 match quality is a separate contract from S21 insertion loss"
    "keep VSWR, return loss, mismatch loss, reflected power, and transmitted power together"
    "use this scalar row as a readable MATLAB optimization objective or constraint"
];

result.status = "ok";
checkNames = fieldnames(checks);
for k = 1:numel(checkNames)
    if ~checks.(checkNames{k})
        result.status = "needs_attention";
        break
    end
end
end
