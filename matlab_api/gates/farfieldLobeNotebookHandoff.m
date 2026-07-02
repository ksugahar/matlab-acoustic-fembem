function result = farfieldLobeNotebookHandoff(thetaValues, phiValues, lobeRow, options)
%farfieldLobeNotebookHandoff Check one far-field lobe row.
%
% A notebook panel should not receive a naked gain value.  The row needs a
% stable lobe identity, frequency, theta/phi location, polarization basis,
% accepted-power normalization, gain/directivity units, and the power identity
% G = eta_rad * D.  This readable helper composes the far-field metadata and
% gain/directivity gates without hiding the assumptions from students.

arguments
    thetaValues (:,:) double
    phiValues (:,:) double
    lobeRow (1,1) struct
    options.FrequencyHz (1,1) double = NaN
    options.ExpectedFrequencyHz (1,1) double = NaN
    options.AngleUnit (1,1) string = "deg"
    options.ExpectedAngleUnit (1,1) string = "deg"
    options.CoordinateSystem (1,1) string = "spherical"
    options.ExpectedCoordinateSystem (1,1) string = "spherical"
    options.PolarizationBasis (1,1) string = "theta_phi"
    options.ExpectedPolarizationBasis (1,1) string = "theta_phi"
    options.Quantity (1,1) string = "gain"
    options.ExpectedQuantity (1,1) string = "gain"
    options.QuantityUnit (1,1) string = "dBi"
    options.ExpectedQuantityUnit (1,1) string = "dBi"
    options.Normalization (1,1) string = "accepted_power"
    options.ExpectedNormalization (1,1) string = "accepted_power"
    options.FieldComponents (:,:) string = ["Etheta" "Ephi"]
    options.RequiredComponents (:,:) string = ["Etheta" "Ephi"]
    options.RequiredPhiValuesDeg (:,:) double = zeros(0, 1)
    options.MinThetaSpanDeg (1,1) double = 180
    options.RowCount (1,1) double = NaN
    options.LobeIdKey (1,1) string = "lobeId"
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

metadata = farfieldPatternMetadata(thetaValues, phiValues, ...
    "FrequencyHz", options.FrequencyHz, ...
    "ExpectedFrequencyHz", options.ExpectedFrequencyHz, ...
    "AngleUnit", options.AngleUnit, ...
    "ExpectedAngleUnit", options.ExpectedAngleUnit, ...
    "CoordinateSystem", options.CoordinateSystem, ...
    "ExpectedCoordinateSystem", options.ExpectedCoordinateSystem, ...
    "PolarizationBasis", options.PolarizationBasis, ...
    "ExpectedPolarizationBasis", options.ExpectedPolarizationBasis, ...
    "Quantity", options.Quantity, ...
    "ExpectedQuantity", options.ExpectedQuantity, ...
    "QuantityUnit", options.QuantityUnit, ...
    "ExpectedQuantityUnit", options.ExpectedQuantityUnit, ...
    "Normalization", options.Normalization, ...
    "ExpectedNormalization", options.ExpectedNormalization, ...
    "FieldComponents", options.FieldComponents, ...
    "RequiredComponents", options.RequiredComponents, ...
    "RequiredPhiValuesDeg", options.RequiredPhiValuesDeg, ...
    "MinThetaSpanDeg", options.MinThetaSpanDeg, ...
    "RowCount", options.RowCount, ...
    "Tolerance", options.Tolerance);

lobeId = getText(lobeRow, [options.LobeIdKey, "lobeId", "lobe_id", "rowId", "row_id", "caseId", "case_id", "label"], "");
rowFrequencyHz = getNumber(lobeRow, ["frequencyHz", "frequency_hz", "frequency_Hz", "freqHz", "freq_hz"], NaN);
thetaDeg = getNumber(lobeRow, ["thetaDeg", "theta_deg", "theta"], NaN);
phiDeg = getNumber(lobeRow, ["phiDeg", "phi_deg", "phi"], NaN);
rowBasis = normalizedBasis(getText(lobeRow, ["polarizationBasis", "polarization_basis", "basis", "polarization"], ""));
rowNormalization = normalizedToken(getText(lobeRow, ["normalization", "powerNormalization", "power_normalization"], ""));
gainUnit = strtrim(getText(lobeRow, ["gainUnit", "gain_unit", "gainQuantityUnit", "gain_quantity_unit"], ""));
directivityUnit = strtrim(getText(lobeRow, ["directivityUnit", "directivity_unit"], ""));
gainDbi = getNumber(lobeRow, ["gainDbi", "gain_dbi", "gain_dBi", "gainDb", "gain_db"], NaN);
directivityDbi = getNumber(lobeRow, ["directivityDbi", "directivity_dbi", "directivity_dBi", "directivityDb", "directivity_db"], NaN);
radiatedPowerW = getNumber(lobeRow, ["radiatedPowerW", "radiated_power_w", "radiated_power_W", "prad_w"], NaN);
acceptedPowerW = getNumber(lobeRow, ["acceptedPowerW", "accepted_power_w", "accepted_power_W", "pacc_w"], NaN);

gainGate = farfieldGainDirectivityEfficiency(gainDbi, directivityDbi, ...
    "RadiatedPowerW", radiatedPowerW, ...
    "AcceptedPowerW", acceptedPowerW, ...
    "Normalization", rowNormalization, ...
    "ExpectedNormalization", options.ExpectedNormalization, ...
    "Tolerance", options.Tolerance);

checks = struct( ...
    "metadataOk", metadata.status == "ok", ...
    "gainDirectivityOk", gainGate.status == "ok", ...
    "lobeIdRecorded", strlength(strtrim(lobeId)) > 0, ...
    "rowFrequencyRecorded", ~isnan(rowFrequencyHz), ...
    "rowFrequencyMatchesMetadata", ~isnan(rowFrequencyHz) && ~isnan(metadata.frequencyHz) && ...
        abs(rowFrequencyHz - metadata.frequencyHz) <= options.Tolerance * max(1, abs(metadata.frequencyHz)), ...
    "thetaRecorded", ~isnan(thetaDeg), ...
    "thetaOnExportGrid", ~isnan(thetaDeg) && containsAngle(metadata.thetaValuesDeg, thetaDeg, options.Tolerance), ...
    "phiRecorded", ~isnan(phiDeg), ...
    "phiOnExportGrid", ~isnan(phiDeg) && containsAngle(metadata.phiValuesDeg, phiDeg, options.Tolerance), ...
    "polarizationBasisRecorded", strlength(rowBasis) > 0, ...
    "polarizationBasisMatchesMetadata", rowBasis == metadata.polarizationBasis, ...
    "normalizationRecorded", strlength(rowNormalization) > 0, ...
    "normalizationMatchesMetadata", rowNormalization == metadata.normalization, ...
    "normalizationIsAcceptedPower", rowNormalization == "accepted_power", ...
    "gainUnitRecorded", strlength(gainUnit) > 0, ...
    "gainUnitIsDbi", gainUnit == "dBi", ...
    "directivityUnitRecorded", strlength(directivityUnit) > 0, ...
    "directivityUnitIsDbi", directivityUnit == "dBi");

result = struct();
result.kind = "farfield_lobe_notebook_handoff";
result.policy = "readable_farfield_lobe_row_before_notebook_panel";
result.lobeId = lobeId;
result.frequencyHz = rowFrequencyHz;
result.thetaDeg = thetaDeg;
result.phiDeg = phiDeg;
result.polarizationBasis = rowBasis;
result.normalization = rowNormalization;
result.gainUnit = gainUnit;
result.directivityUnit = directivityUnit;
result.metadata = metadata;
result.gainDirectivity = gainGate;
result.checks = checks;
result.notes = [
    "run this after farfieldPatternMetadata and before plotting or ranking lobes"
    "a notebook row should expose lobe identity, angular location, polarization, normalization, and G = eta_rad * D together"
    "negative controls should include missing lobe id, phi outside the cut grid, and gain above directivity"
];

result.status = "ok";
checkValues = struct2cell(checks);
if ~all([checkValues{:}])
    result.status = "needs_attention";
end
end


function value = getText(row, names, defaultValue)
value = defaultValue;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(row, name)
        raw = row.(name);
        if isstring(raw) || ischar(raw)
            value = string(raw);
        elseif isnumeric(raw) || islogical(raw)
            value = string(raw);
        end
        return
    end
end
value = string(value);
end


function value = getNumber(row, names, defaultValue)
value = defaultValue;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(row, name)
        raw = row.(name);
        if isnumeric(raw) || islogical(raw)
            value = double(raw);
        elseif isstring(raw) || ischar(raw)
            value = str2double(raw);
        end
        return
    end
end
end


function ok = containsAngle(values, target, tolerance)
ok = any(abs(values(:) - target) <= max(tolerance, 1e-9));
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
