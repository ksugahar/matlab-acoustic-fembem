function result = touchstonePortMetadata(ports, options)
%touchstonePortMetadata Check Touchstone table metadata before rows.
%
% A passive S-parameter row is not enough if the table has lost its port
% names, port order, RI/MA/DB format, frequency unit, network parameter, or
% reference impedance.  This helper keeps those small contracts visible before
% a row is reused in optimization, equivalent-circuit fitting, or FEM/BEM
% notebooks.

arguments
    ports (:,:) string
    options.RequiredPorts (:,:) string = ["P1" "P2"]
    options.PortOrder (:,:) string = strings(0, 1)
    options.NetworkParameter (1,1) string = ""
    options.ExpectedNetworkParameter (1,1) string = "S"
    options.DataFormat (1,1) string = ""
    options.ExpectedDataFormat (1,1) string = ""
    options.FrequencyUnit (1,1) string = ""
    options.ExpectedFrequencyUnit (1,1) string = ""
    options.Z0 (1,1) double = NaN
    options.ExpectedZ0 (1,1) double = NaN
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-12
end

portNames = strip(ports(:));
portNames(portNames == "") = [];
requiredPorts = strip(options.RequiredPorts(:));
requiredPorts(requiredPorts == "") = [];
expectedOrder = strip(options.PortOrder(:));
expectedOrder(expectedOrder == "") = [];

network = upper(strtrim(options.NetworkParameter));
expectedNetwork = upper(strtrim(options.ExpectedNetworkParameter));
dataFormat = upper(strtrim(options.DataFormat));
expectedFormat = upper(strtrim(options.ExpectedDataFormat));
frequencyUnit = normalizedFrequencyUnit(options.FrequencyUnit);
expectedFrequencyUnit = normalizedFrequencyUnit(options.ExpectedFrequencyUnit);
z0 = options.Z0;
expectedZ0 = options.ExpectedZ0;

portsUnique = numel(unique(portNames)) == numel(portNames);
requiredPresent = all(ismember(requiredPorts, portNames));
if isempty(expectedOrder)
    orderMatches = true;
else
    orderMatches = numel(portNames) >= numel(expectedOrder) && ...
        all(portNames(1:numel(expectedOrder)) == expectedOrder);
end

networkAllowed = any(network == ["S", "Y", "Z", "H", "G"]);
if strlength(expectedNetwork) > 0
    networkMatches = networkAllowed && network == expectedNetwork;
else
    networkMatches = networkAllowed;
end

formatAllowed = any(dataFormat == ["RI", "MA", "DB"]);
if strlength(expectedFormat) > 0
    formatMatches = formatAllowed && dataFormat == expectedFormat;
else
    formatMatches = formatAllowed;
end

frequencyAllowed = strlength(frequencyUnit) > 0;
if strlength(expectedFrequencyUnit) > 0
    frequencyMatches = frequencyAllowed && frequencyUnit == expectedFrequencyUnit;
else
    frequencyMatches = frequencyAllowed;
end

z0Recorded = ~isnan(z0);
z0Positive = z0Recorded && z0 > 0;
if isnan(expectedZ0)
    z0Matches = z0Positive;
else
    z0Matches = z0Positive && abs(z0 - expectedZ0) <= options.Tolerance * max(1, abs(expectedZ0));
end

checks = struct( ...
    "portNamesRecorded", ~isempty(portNames), ...
    "portsUnique", portsUnique, ...
    "requiredPortsPresent", requiredPresent, ...
    "portOrderMatchesExpected", orderMatches, ...
    "networkParameterRecorded", strlength(network) > 0, ...
    "networkParameterMatchesExpected", networkMatches, ...
    "touchstoneFormatRecorded", strlength(dataFormat) > 0, ...
    "touchstoneFormatMatchesExpected", formatMatches, ...
    "frequencyUnitRecorded", strlength(frequencyUnit) > 0, ...
    "frequencyUnitMatchesExpected", frequencyMatches, ...
    "referenceImpedanceRecorded", z0Recorded, ...
    "referenceImpedanceMatchesExpected", z0Matches);

result = struct();
result.kind = "touchstone_port_metadata";
result.policy = "readable_touchstone_metadata_before_row_values";
result.ports = portNames;
result.portCount = numel(portNames);
result.requiredPorts = requiredPorts;
result.expectedPortOrder = expectedOrder;
result.networkParameter = network;
result.dataFormat = dataFormat;
result.frequencyUnit = frequencyUnit;
result.z0 = z0;
result.checks = checks;
result.notes = [
    "check metadata before RI/MA/DB row values are converted to complex S-parameters"
    "a passive and reciprocal row can still be wrong if ports are swapped or z0 is missing"
    "this is the MATLAB teaching companion to CST/radia Touchstone port metadata gates"
];

result.status = "ok";
checkValues = struct2cell(checks);
if ~all([checkValues{:}])
    result.status = "needs_attention";
end
end


function unit = normalizedFrequencyUnit(value)
key = lower(strtrim(value));
switch key
    case "hz"
        unit = "Hz";
    case "khz"
        unit = "kHz";
    case "mhz"
        unit = "MHz";
    case "ghz"
        unit = "GHz";
    otherwise
        unit = "";
end
end
