function result = farfieldPatternMetadata(thetaValues, phiValues, options)
%farfieldPatternMetadata Check far-field table metadata first.
%
% Far-field gain, directivity, RCS, or Etheta/Ephi rows are not evidence until
% the table says which frequency, angle unit, theta/phi cut, polarization
% basis, quantity unit, row count, and power normalization produced them.  This
% helper is intentionally small so the metadata contract stays visible in
% MATLAB notebooks before a Gypsilab/NGSolve.BEM comparison or optimization.

arguments
    thetaValues (:,:) double
    phiValues (:,:) double
    options.FrequencyHz (1,1) double = NaN
    options.ExpectedFrequencyHz (1,1) double = NaN
    options.AngleUnit (1,1) string = ""
    options.ExpectedAngleUnit (1,1) string = "deg"
    options.CoordinateSystem (1,1) string = ""
    options.ExpectedCoordinateSystem (1,1) string = "spherical"
    options.PolarizationBasis (1,1) string = ""
    options.ExpectedPolarizationBasis (1,1) string = "theta_phi"
    options.Quantity (1,1) string = ""
    options.ExpectedQuantity (1,1) string = "gain"
    options.QuantityUnit (1,1) string = ""
    options.ExpectedQuantityUnit (1,1) string = "dBi"
    options.Normalization (1,1) string = ""
    options.ExpectedNormalization (1,1) string = "accepted_power"
    options.FieldComponents (:,:) string = strings(0, 1)
    options.RequiredComponents (:,:) string = ["Etheta" "Ephi"]
    options.RequiredPhiValuesDeg (:,:) double = zeros(0, 1)
    options.MinThetaSpanDeg (1,1) double = 180
    options.RowCount (1,1) double = NaN
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

angleUnit = normalizedAngleUnit(options.AngleUnit);
expectedAngleUnit = normalizedAngleUnit(options.ExpectedAngleUnit);
thetaDeg = convertAngles(thetaValues(:), angleUnit);
phiDeg = convertAngles(phiValues(:), angleUnit);

components = strip(options.FieldComponents(:));
components(components == "") = [];
required = strip(options.RequiredComponents(:));
required(required == "") = [];

coordinateSystem = normalizedToken(options.CoordinateSystem);
expectedCoordinateSystem = normalizedToken(options.ExpectedCoordinateSystem);
polarizationBasis = normalizedBasis(options.PolarizationBasis);
expectedPolarizationBasis = normalizedBasis(options.ExpectedPolarizationBasis);
quantity = normalizedToken(options.Quantity);
expectedQuantity = normalizedToken(options.ExpectedQuantity);
normalization = normalizedToken(options.Normalization);
expectedNormalization = normalizedToken(options.ExpectedNormalization);
quantityUnit = strtrim(options.QuantityUnit);
expectedQuantityUnit = strtrim(options.ExpectedQuantityUnit);

frequencyRecorded = ~isnan(options.FrequencyHz);
frequencyPositive = frequencyRecorded && options.FrequencyHz > 0;
if isnan(options.ExpectedFrequencyHz)
    frequencyMatches = true;
else
    frequencyMatches = frequencyRecorded && ...
        abs(options.FrequencyHz - options.ExpectedFrequencyHz) <= ...
        options.Tolerance * max(1, abs(options.ExpectedFrequencyHz));
end

thetaSpanDeg = NaN;
thetaRangeOk = false;
thetaStrict = true;
if ~isempty(thetaDeg)
    thetaSpanDeg = max(thetaDeg) - min(thetaDeg);
    thetaRangeOk = all(isfinite(thetaDeg)) && min(thetaDeg) >= -options.Tolerance && ...
        max(thetaDeg) <= 180 + options.Tolerance;
    thetaStrict = all(diff(thetaDeg) > 0);
end

phiRangeOk = false;
phiUnique = true;
if ~isempty(phiDeg)
    phiRangeOk = all(isfinite(phiDeg)) && min(phiDeg) >= -360 - options.Tolerance && ...
        max(phiDeg) <= 360 + options.Tolerance;
    phiUnique = numel(unique(round(phiDeg, 12))) == numel(phiDeg);
end

requiredPhi = options.RequiredPhiValuesDeg(:);
requiredPhiPresent = true;
for k = 1:numel(requiredPhi)
    requiredPhiPresent = requiredPhiPresent && any(abs(phiDeg - requiredPhi(k)) <= max(options.Tolerance, 1e-9));
end

rowCount = options.RowCount;
expectedGridRows = numel(thetaDeg) * max(1, numel(phiDeg));

checks = struct( ...
    "frequencyRecorded", frequencyRecorded, ...
    "frequencyPositive", frequencyPositive, ...
    "frequencyMatchesExpected", frequencyMatches, ...
    "angleUnitRecorded", strlength(angleUnit) > 0, ...
    "angleUnitMatchesExpected", angleUnit == expectedAngleUnit, ...
    "coordinateSystemRecorded", strlength(coordinateSystem) > 0, ...
    "coordinateSystemMatchesExpected", coordinateSystem == expectedCoordinateSystem, ...
    "polarizationBasisRecorded", strlength(polarizationBasis) > 0, ...
    "polarizationBasisMatchesExpected", polarizationBasis == expectedPolarizationBasis, ...
    "quantityRecorded", strlength(quantity) > 0, ...
    "quantityMatchesExpected", quantity == expectedQuantity, ...
    "quantityUnitRecorded", strlength(quantityUnit) > 0, ...
    "quantityUnitMatchesExpected", quantityUnit == expectedQuantityUnit, ...
    "normalizationRecorded", strlength(normalization) > 0, ...
    "normalizationMatchesExpected", normalization == expectedNormalization, ...
    "fieldComponentsRecorded", ~isempty(components), ...
    "requiredComponentsPresent", all(ismember(required, components)), ...
    "thetaGridRecorded", ~isempty(thetaDeg), ...
    "thetaGridStrictlyIncreasing", thetaStrict, ...
    "thetaRangeDegreesOk", thetaRangeOk, ...
    "thetaSpanOk", ~isnan(thetaSpanDeg) && thetaSpanDeg + options.Tolerance >= options.MinThetaSpanDeg, ...
    "phiGridRecorded", ~isempty(phiDeg), ...
    "phiValuesUnique", phiUnique, ...
    "phiRangeDegreesOk", phiRangeOk, ...
    "requiredPhiValuesPresent", requiredPhiPresent, ...
    "rowCountRecorded", ~isnan(rowCount), ...
    "rowCountCoversGrid", ~isnan(rowCount) && rowCount >= expectedGridRows);

result = struct();
result.kind = "farfield_pattern_metadata";
result.policy = "readable_farfield_metadata_before_pattern_values";
result.frequencyHz = options.FrequencyHz;
result.expectedFrequencyHz = options.ExpectedFrequencyHz;
result.angleUnit = angleUnit;
result.thetaValuesDeg = thetaDeg;
result.phiValuesDeg = phiDeg;
result.thetaSpanDeg = thetaSpanDeg;
result.coordinateSystem = coordinateSystem;
result.polarizationBasis = polarizationBasis;
result.quantity = quantity;
result.quantityUnit = quantityUnit;
result.normalization = normalization;
result.fieldComponents = components;
result.requiredComponents = required;
result.requiredPhiValuesDeg = requiredPhi;
result.rowCount = rowCount;
result.expectedGridRows = expectedGridRows;
result.checks = checks;
result.notes = [
    "metadata comes before far-field gain, directivity, RCS, or Etheta/Ephi row values"
    "angle units, theta/phi cuts, polarization basis, and power normalization are part of the physics contract"
    "this is the MATLAB teaching companion to CST/radia far-field metadata gates"
];

result.status = "ok";
checkValues = struct2cell(checks);
if ~all([checkValues{:}])
    result.status = "needs_attention";
end
end


function unit = normalizedAngleUnit(value)
key = lower(strtrim(value));
switch key
    case {"deg", "degree", "degrees"}
        unit = "deg";
    case {"rad", "radian", "radians"}
        unit = "rad";
    otherwise
        unit = "";
end
end


function valuesDeg = convertAngles(values, unit)
if unit == "rad"
    valuesDeg = rad2deg(values);
else
    valuesDeg = values;
end
end


function token = normalizedToken(value)
token = lower(strrep(strrep(strtrim(value), "-", "_"), "/", "_"));
end


function token = normalizedBasis(value)
raw = normalizedToken(value);
switch raw
    case {"theta_phi", "etheta_ephi", "spherical_theta_phi"}
        token = "theta_phi";
    case {"lhcp_rhcp", "rhcp_lhcp"}
        token = "circular";
    otherwise
        token = raw;
end
end
