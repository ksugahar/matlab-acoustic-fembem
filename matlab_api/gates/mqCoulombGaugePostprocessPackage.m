function result = mqCoulombGaugePostprocessPackage(artifacts, options)
%mqCoulombGaugePostprocessPackage Bundle MQS E-field postprocess gates.
%
% This readable helper keeps a magneto-quasistatic A-phi solve, Coulomb-gauge
% postprocess, spatial-potential boundary condition, recovered E-field row, and
% Darwin/full-wave validity envelope together before a notebook presents the
% electric field as evidence.

arguments
    artifacts (1,:) struct
    options.ExpectedCaseId (1,1) string = ""
    options.ExpectedMeshId (1,1) string = ""
    options.ExpectedFrequencyHz (1,1) double = NaN
    options.MaxFrequencyRatioToFullwave (1,1) double {mustBeNonnegative} = 0.1
    options.FrequencyTolerance (1,1) double {mustBeNonnegative} = 1e-12
end

requiredKinds = [
    "mqs_solution"
    "coulomb_gauge"
    "spatial_potential"
    "electric_field"
    "validity_envelope"
];

n = numel(artifacts);
kinds = strings(n, 1);
caseIds = strings(n, 1);
meshIds = strings(n, 1);
paths = strings(n, 1);
statuses = strings(n, 1);
policies = strings(n, 1);
frequencies = NaN(n, 1);

for k = 1:n
    row = artifacts(k);
    kinds(k) = normalizeToken(stringField(row, ["kind", "artifactKind", "type"]));
    caseIds(k) = stringField(row, ["caseId", "case_id", "modelId"]);
    meshIds(k) = stringField(row, ["meshId", "mesh_id", "volId"]);
    paths(k) = stringField(row, ["path", "file", "artifactPath"]);
    statuses(k) = normalizeToken(stringField(row, ["status", "gateStatus", "validationStatus"]));
    policies(k) = normalizeToken(stringField(row, ["gatePolicy", "gate_policy", "policy", "validator"]));
    frequencies(k) = doubleField(row, ["frequencyHz", "frequency_Hz", "frequency_hz"]);
end

presentKinds = unique(kinds(kinds ~= ""));
frequencyValues = frequencies(~isnan(frequencies));
frequencyRelSpan = 0;
if ~isempty(frequencyValues)
    frequencyRelSpan = (max(frequencyValues) - min(frequencyValues)) / max(max(abs(frequencyValues)), 1);
end

mqsRows = kinds == "mqs_solution";
gaugeRows = kinds == "coulomb_gauge";
potentialRows = kinds == "spatial_potential";
efieldRows = kinds == "electric_field";
validityRows = kinds == "validity_envelope";

checks = struct( ...
    "requiredKindsPresent", all(ismember(requiredKinds, presentKinds)), ...
    "caseIdsPresent", all(caseIds ~= ""), ...
    "caseIdsUnique", isscalar(unique(caseIds(caseIds ~= ""))), ...
    "meshIdsPresent", all(meshIds ~= ""), ...
    "meshIdsUnique", isscalar(unique(meshIds(meshIds ~= ""))), ...
    "frequenciesPresent", all(~isnan(frequencies)), ...
    "frequenciesMatch", ~isempty(frequencyValues) && frequencyRelSpan <= options.FrequencyTolerance, ...
    "pathsPresent", all(paths ~= ""), ...
    "upstreamStatusOk", all(ismember(statuses, ["ok", "pass", "passed", "verified"])), ...
    "gatePoliciesKnown", policiesKnown(kinds, policies), ...
    "mqsAPhiFormulationRecorded", any(mqsRows & ismember(arrayfun(@(x) normalizeToken(stringField(x, "formulation")), artifacts).', ["a_phi", "a_phi_mqs", "mqs_a_phi"])), ...
    "coulombGaugeConditionRecorded", any(gaugeRows & ismember(arrayfun(@(x) normalizeToken(stringField(x, "gaugeCondition")), artifacts).', ["coulomb", "div_a_zero", "divergence_a_zero"])), ...
    "spatialPotentialBcRecorded", any(potentialRows & arrayfun(@(x) stringField(x, "boundaryConditionSource") ~= "", artifacts).'), ...
    "electricFieldUnitRecorded", any(efieldRows & ismember(arrayfun(@(x) normalizeToken(stringField(x, ["EUnit", "eUnit"])), artifacts).', ["v_per_m", "v/m"])), ...
    "validityEnvelopeOk", validityEnvelopeOk(artifacts(validityRows), options.MaxFrequencyRatioToFullwave));

if options.ExpectedCaseId ~= ""
    checks.expectedCaseIdMatches = all(caseIds(caseIds ~= "") == options.ExpectedCaseId);
end
if options.ExpectedMeshId ~= ""
    checks.expectedMeshIdMatches = all(meshIds(meshIds ~= "") == options.ExpectedMeshId);
end
if ~isnan(options.ExpectedFrequencyHz)
    checks.expectedFrequencyMatches = ~isempty(frequencyValues) && ...
        max(abs(frequencyValues - options.ExpectedFrequencyHz) ./ max(abs(frequencyValues), abs(options.ExpectedFrequencyHz))) <= options.FrequencyTolerance;
end

result = struct();
result.kind = "mqs_coulomb_gauge_efield_postprocess_package";
result.policy = "readable_mqs_coulomb_gauge_efield_postprocess_gate";
result.requiredKinds = requiredKinds;
result.presentKinds = presentKinds;
result.caseIds = unique(caseIds(caseIds ~= ""));
result.meshIds = unique(meshIds(meshIds ~= ""));
result.frequenciesHz = unique(frequencyValues);
result.maxFrequencyRelSpan = frequencyRelSpan;
result.maxFrequencyRatioToFullwave = options.MaxFrequencyRatioToFullwave;
result.checks = checks;
result.notes = [
    "A recovered E-field from MQS is a postprocess package, not an isolated field row."
    "Carry the A-phi formulation, Coulomb gauge condition, conductor-surface potential boundary source, E-field units, and Darwin/full-wave validity envelope together."
    "This is a readable MATLAB companion to radia-mcp mqs_coulomb_gauge_efield_postprocess_gate."
];

result.status = "ok";
names = fieldnames(checks);
for k = 1:numel(names)
    if ~checks.(names{k})
        result.status = "needs_attention";
        break
    end
end
end


function ok = policiesKnown(kinds, policies)
ok = true;
for k = 1:numel(kinds)
    switch kinds(k)
        case "mqs_solution"
            allowed = ["mqs_a_phi_solution", "a_phi_mqs_solution", "mqs_solution"];
        case "coulomb_gauge"
            allowed = ["coulomb_gauge_postprocess", "coulomb_gauge_condition"];
        case "spatial_potential"
            allowed = ["electrostatic_potential_postprocess", "spatial_potential_solve"];
        case "electric_field"
            allowed = ["efield_gradient_recovery", "electric_field_postprocess"];
        case "validity_envelope"
            allowed = ["mqs_darwin_fullwave_validity_envelope", "quasistatic_validity_envelope"];
        otherwise
            ok = false;
            return
    end
    if ~ismember(policies(k), allowed)
        ok = false;
        return
    end
end
end


function ok = validityEnvelopeOk(rows, maxRatio)
ok = false;
if isempty(rows)
    return
end
for k = 1:numel(rows)
    ratio = doubleField(rows(k), "frequencyRatioToFullwaveLimit");
    dominant = logicalField(rows(k), "dominantInductive");
    reference = normalizeToken(stringField(rows(k), "comparisonReference"));
    ok = ~isnan(ratio) && ratio <= maxRatio && dominant && ...
        ismember(reference, ["darwin", "full_wave", "darwin_full_wave", "fullwave"]);
    if ok
        return
    end
end
end


function value = stringField(row, names)
names = string(names);
value = "";
for k = 1:numel(names)
    name = char(names(k));
    if isfield(row, name)
        raw = row.(name);
        if ~(isstring(raw) && isscalar(raw) && raw == "")
            value = string(raw);
            return
        end
    end
end
end


function value = doubleField(row, names)
raw = stringField(row, names);
if raw == ""
    value = NaN;
else
    value = str2double(raw);
end
end


function value = logicalField(row, names)
raw = stringField(row, names);
switch normalizeToken(raw)
    case {"true", "1", "yes"}
        value = true;
    otherwise
        value = false;
end
end


function token = normalizeToken(value)
token = lower(strtrim(string(value)));
token = replace(token, "-", "_");
token = replace(token, " ", "_");
end
