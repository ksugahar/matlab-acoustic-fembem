function result = farfieldGainDirectivityEfficiency(gainDbi, directivityDbi, options)
%farfieldGainDirectivityEfficiency Check G = eta_rad * D.
%
% After far-field metadata confirms accepted-power normalization, a gain row
% should still agree with directivity and radiation efficiency:
%
%   eta_rad = P_rad / P_acc,    G = eta_rad * D.
%
% This small readable gate is the MATLAB teaching companion to the CST source
% lane and radia-mcp antenna helper.  It intentionally keeps accepted power,
% radiated power, gain, and directivity in one struct before rows are used by
% notebooks, optimization, or FEM/BEM comparison.

arguments
    gainDbi (1,1) double
    directivityDbi (1,1) double
    options.RadiatedPowerW (1,1) double = NaN
    options.AcceptedPowerW (1,1) double = NaN
    options.Normalization (1,1) string = ""
    options.ExpectedNormalization (1,1) string = "accepted_power"
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

gainLinear = 10^(gainDbi / 10);
directivityLinear = 10^(directivityDbi / 10);
normalization = normalizedToken(options.Normalization);
expectedNormalization = normalizedToken(options.ExpectedNormalization);

radiatedRecorded = ~isnan(options.RadiatedPowerW);
acceptedRecorded = ~isnan(options.AcceptedPowerW);
radiatedNonnegative = radiatedRecorded && options.RadiatedPowerW >= 0;
acceptedPositive = acceptedRecorded && options.AcceptedPowerW > 0;

radiationEfficiency = NaN;
expectedGainLinear = NaN;
gainRelativeError = NaN;
efficiencyFromGain = gainLinear / directivityLinear;
if acceptedPositive
    radiationEfficiency = options.RadiatedPowerW / options.AcceptedPowerW;
    expectedGainLinear = directivityLinear * radiationEfficiency;
    gainRelativeError = abs(gainLinear - expectedGainLinear) / ...
        max([abs(expectedGainLinear), abs(gainLinear), 1]);
end

checks = struct( ...
    "normalizationRecorded", strlength(normalization) > 0, ...
    "normalizationMatchesExpected", normalization == expectedNormalization, ...
    "radiatedPowerRecorded", radiatedRecorded, ...
    "acceptedPowerRecorded", acceptedRecorded, ...
    "radiatedPowerNonnegative", radiatedNonnegative, ...
    "acceptedPowerPositive", acceptedPositive, ...
    "radiationEfficiencyInZeroOne", ~isnan(radiationEfficiency) && ...
        radiationEfficiency >= 0 && radiationEfficiency <= 1 + options.Tolerance, ...
    "gainNotAboveDirectivity", gainLinear <= directivityLinear + options.Tolerance, ...
    "gainMatchesDirectivityTimesEfficiency", ~isnan(gainRelativeError) && ...
        gainRelativeError <= options.Tolerance);

result = struct();
result.kind = "farfield_gain_directivity_efficiency";
result.policy = "readable_farfield_gain_after_accepted_power_metadata";
result.gainDbi = gainDbi;
result.directivityDbi = directivityDbi;
result.gainLinear = gainLinear;
result.directivityLinear = directivityLinear;
result.radiatedPowerW = options.RadiatedPowerW;
result.acceptedPowerW = options.AcceptedPowerW;
result.normalization = normalization;
result.radiationEfficiency = radiationEfficiency;
result.efficiencyFromGainOverDirectivity = efficiencyFromGain;
result.expectedGainLinear = expectedGainLinear;
result.gainRelativeError = gainRelativeError;
result.checks = checks;
result.notes = [
    "run farfieldPatternMetadata before this numeric row gate"
    "accepted-power gain should equal radiation efficiency times directivity"
    "gain above directivity usually means normalization, unit, or row selection is wrong"
];

result.status = "ok";
checkValues = struct2cell(checks);
if ~all([checkValues{:}])
    result.status = "needs_attention";
end
end


function token = normalizedToken(value)
token = lower(strrep(strrep(strtrim(value), "-", "_"), "/", "_"));
end
