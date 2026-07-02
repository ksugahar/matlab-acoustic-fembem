function result = educationalTouchstoneSolverReadyPreflight(s11, s21, options)
%EDUCATIONALTOUCHSTONESOLVERREADYPREFLIGHT Bundle Touchstone row gates.
%
% The inputs are already-normalized complex S-parameters.  This helper keeps
% the row-level metadata and algebraic checks together before a CST/Touchstone,
% ngsolve.bem, or measurement row is reused as a MATLAB optimization objective,
% constraint, or notebook result.

arguments
    s11 (1,1) double
    s21 (1,1) double
    options.S12 (1,1) double = NaN
    options.S22 (1,1) double = NaN
    options.Frequency (1,1) double = NaN
    options.FrequencyUnit (1,1) string = "GHz"
    options.DataFormat (1,1) string = "MA"
    options.Z0 (1,1) double {mustBePositive} = 50
    options.ReturnLossMinDb (1,1) double = NaN
    options.VswrMax (1,1) double = NaN
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

s12 = options.S12;
s22 = options.S22;
if isnan(s12)
    s12 = s21;
end
if isnan(s22)
    s22 = s11;
end

frequencyRecorded = ~isnan(options.Frequency);
dataFormat = upper(strtrim(options.DataFormat));
formatRecorded = any(dataFormat == ["RI", "MA", "DB"]);
frequencyScale = frequencyUnitScale(options.FrequencyUnit);
frequencyHz = NaN;
if frequencyRecorded
    frequencyHz = options.Frequency * frequencyScale;
end

equivalent = educationalTouchstoneEquivalentCircuit(s11, s21, ...
    "S12", s12, ...
    "S22", s22, ...
    "Z0", options.Z0, ...
    "Tolerance", options.Tolerance);
match = educationalTouchstonePortMatch(s11, ...
    "ReturnLossMinDb", options.ReturnLossMinDb, ...
    "VswrMax", options.VswrMax, ...
    "Tolerance", options.Tolerance);

checks = struct( ...
    "frequencyRecorded", frequencyRecorded, ...
    "formatRecorded", formatRecorded, ...
    "referenceImpedanceContractOk", equivalent.status == "ok", ...
    "sparameterPassivityOk", equivalent.checks.sparameterPassivityOk, ...
    "sparameterReciprocityOk", equivalent.checks.sparameterReciprocityOk, ...
    "portMatchContractOk", match.status == "ok", ...
    "returnLossLimitOk", match.checks.returnLossLimitOk, ...
    "vswrLimitOk", match.checks.vswrLimitOk);

result = struct();
result.kind = "educational_touchstone_solver_ready_preflight";
result.policy = "readable_touchstone_solver_ready_row_preflight";
result.frequency = options.Frequency;
result.frequencyUnit = options.FrequencyUnit;
result.frequencyHz = frequencyHz;
result.dataFormat = dataFormat;
result.z0 = options.Z0;
result.s = struct("s11", s11, "s21", s21, "s12", s12, "s22", s22);
result.equivalent = equivalent;
result.portMatch = match;
result.checks = checks;
result.notes = [
    "normalize RI/MA/DB to complex S-parameters before calling this helper"
    "do not reuse a row in optimization unless frequency, format, Z0, passivity, reciprocity, and S11 match are carried together"
    "this is a readable MATLAB companion to CST/radia Touchstone solver-ready preflight gates"
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


function scale = frequencyUnitScale(unit)
switch lower(strtrim(unit))
    case "hz"
        scale = 1;
    case "khz"
        scale = 1e3;
    case "mhz"
        scale = 1e6;
    case "ghz"
        scale = 1e9;
    otherwise
        error("educationalTouchstoneSolverReadyPreflight:frequencyUnit", ...
            "Unsupported frequency unit: %s", unit);
end
end
