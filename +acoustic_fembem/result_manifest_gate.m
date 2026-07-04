function report = result_manifest_gate(result, options)
%RESULT_MANIFEST_GATE Check notebook-ready result JSON/struct metadata.
%
% report = acoustic_fembem.result_manifest_gate(result)
% report = acoustic_fembem.result_manifest_gate(jsonFile)
%
% A reusable CAE result should carry enough provenance for a notebook panel:
% run date, MATLAB/Radia versions, pass/fail status, and a compact timing
% breakdown.  Keep timing to the heaviest few stages so students can scan it.

arguments
    result
    options.RequiredVersions (1,:) string = ["matlab"]
    options.RequireRadiaVersion (1,1) logical = true
    options.RequiredTimingStages (1,:) string = strings(1, 0)
    options.MinTimingStages (1,1) double = 1
    options.MaxTimingStages (1,1) double = 4
    options.RequireRunDate (1,1) logical = true
    options.ExpectedCreatedAtUtc (1,1) string = ""
    options.ExpectedRunDateUtc (1,1) string = ""
    options.MaxCreatedRunSkewSeconds (1,1) double = Inf
    options.RequireExecutionSessionId (1,1) logical = false
    options.ExpectedExecutionSessionId (1,1) string = ""
    options.RequireParameterSetArtifact (1,1) logical = false
    options.ExpectedParameterSetArtifactId (1,1) string = ""
    options.ExpectedParameterSetDigest (1,1) string = ""
    options.ExpectedParameterSetPath (1,1) string = ""
    options.ExpectedObjectiveObservableId (1,1) string = ""
    options.ExpectedObjectiveObservableFamily (1,1) string = ""
    options.RequireResultOutputSchema (1,1) logical = false
    options.ExpectedResultOutputSchemaId (1,1) string = ""
    options.ExpectedResultOutputColumns (1,:) string = strings(1, 0)
    options.ExpectedResultOutputUnits = struct()
    options.RequirePhysicsConventionSchema (1,1) logical = false
    options.ExpectedPhysicsConventionSchemaId (1,1) string = ""
    options.RequirePostprocessRowConventionSchema (1,1) logical = false
    options.ExpectedPostprocessRowConventionSchemaId (1,1) string = ""
end

artifact = loadArtifact(result);
if options.MaxCreatedRunSkewSeconds < 0
    error("result_manifest_gate:InvalidOption", "MaxCreatedRunSkewSeconds must be nonnegative")
end
versions = getStructField(artifact, "versions", struct());
timing = getStructField(artifact, "timing_breakdown_s", struct());
execution = getStructField(artifact, "execution", struct());
optimization = getStructField(artifact, "optimization", struct());
resultBlock = getStructField(artifact, "result", struct());
outputBlock = getStructField(artifact, "output", getStructField(artifact, "outputs", struct()));
if ~isstruct(resultBlock)
    resultBlock = struct();
end
if ~isstruct(outputBlock)
    outputBlock = struct();
end

createdAt = string(getStructField(artifact, "created_at_utc", ""));
runDate = string(getStructField(execution, "run_date_utc", getStructField(artifact, "run_date_utc", createdAt)));
runDateRecorded = hasNonemptyField(execution, "run_date_utc") || hasNonemptyField(artifact, "run_date_utc");
createdSeconds = isoDateTimeSeconds(createdAt);
runSeconds = isoDateTimeSeconds(runDate);
createdRunSkewSeconds = NaN;
if ~isnan(createdSeconds) && ~isnan(runSeconds)
    createdRunSkewSeconds = abs(runSeconds - createdSeconds);
end

versionKeys = string(fieldnames(versions));
missingVersions = strings(1, 0);
for key = options.RequiredVersions
    if ~hasNonemptyField(versions, key)
        missingVersions(end + 1) = key; %#ok<AGROW>
    end
end

radiaVersionRecorded = hasAnyNonemptyField(versions, ["radia", "radia_mcp", "radia-ngsolve", "radia_ngsolve"]) || ...
    hasAnyNonemptyField(artifact, ["radia_version", "radia_mcp_version", "radia_ngsolve_version"]);

[timingValues, timingValuesOk] = timingStructToTable(timing);
timingStages = string(timingValues.stage);
missingTimingStages = strings(1, 0);
for key = options.RequiredTimingStages
    if ~any(timingStages == key)
        missingTimingStages(end + 1) = key; %#ok<AGROW>
    end
end

timingValues = sortrows(timingValues, "seconds", "descend");
dominantCount = min(height(timingValues), 4);
dominantTiming = timingValues(1:dominantCount, :);

parameterSetArtifactId = firstStringField( ...
    {artifact, execution, optimization}, ...
    ["parameter_set_artifact_id", "parameterSetArtifactId", ...
     "initial_value_artifact_id", "initialValueArtifactId"]);
parameterSetArtifactIdValues = uniqueStringFields( ...
    {artifact, execution, optimization}, ...
    ["parameter_set_artifact_id", "parameterSetArtifactId", ...
     "initial_value_artifact_id", "initialValueArtifactId"]);
parameterSetDigest = firstStringField( ...
    {artifact, execution, optimization}, ...
    ["parameter_set_digest", "parameterSetDigest", ...
     "parameter_set_sha256", "parameterSetSha256", ...
     "initial_value_digest", "initialValueDigest"]);
parameterSetDigestValues = uniqueStringFields( ...
    {artifact, execution, optimization}, ...
    ["parameter_set_digest", "parameterSetDigest", ...
     "parameter_set_sha256", "parameterSetSha256", ...
     "initial_value_digest", "initialValueDigest"]);
parameterSetPath = firstStringField( ...
    {artifact, execution, optimization}, ...
    ["parameter_set_path", "parameterSetPath", ...
     "initial_value_path", "initialValuePath"]);
parameterSetPathValues = uniqueStringFields( ...
    {artifact, execution, optimization}, ...
    ["parameter_set_path", "parameterSetPath", ...
     "initial_value_path", "initialValuePath"]);
objectiveObservableId = firstStringField( ...
    {artifact, execution, optimization}, ...
    ["objective_observable_id", "objectiveObservableId", ...
     "objective_artifact_id", "objectiveArtifactId"]);
objectiveObservableIdValues = uniqueStringFields( ...
    {artifact, execution, optimization}, ...
    ["objective_observable_id", "objectiveObservableId", ...
     "objective_artifact_id", "objectiveArtifactId"]);
objectiveObservableFamily = firstStringField( ...
    {artifact, execution, optimization}, ...
    ["objective_observable_family", "objectiveObservableFamily", ...
     "objective_family", "objectiveFamily"]);
objectiveObservableFamilyValues = uniqueStringFields( ...
    {artifact, execution, optimization}, ...
    ["objective_observable_family", "objectiveObservableFamily", ...
     "objective_family", "objectiveFamily"]);
executionSessionId = firstStringField( ...
    {execution, artifact}, ...
    ["execution_session_id", "executionSessionId", ...
     "session_id", "sessionId", ...
     "shared_engine", "sharedEngine", ...
     "matlab_engine_session", "matlabEngineSession", ...
     "engine_session", "engineSession"]);
executionSessionIdValues = uniqueStringFields( ...
    {execution, artifact}, ...
    ["execution_session_id", "executionSessionId", ...
     "session_id", "sessionId", ...
     "shared_engine", "sharedEngine", ...
     "matlab_engine_session", "matlabEngineSession", ...
     "engine_session", "engineSession"]);
resultOutputSchemaId = firstStringField( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["result_output_schema_id", "resultOutputSchemaId", ...
     "output_schema_id", "outputSchemaId", ...
     "table_schema_id", "tableSchemaId"]);
resultOutputSchemaIdValues = uniqueStringFields( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["result_output_schema_id", "resultOutputSchemaId", ...
     "output_schema_id", "outputSchemaId", ...
     "table_schema_id", "tableSchemaId"]);
resultOutputColumns = firstStringListField( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["result_output_columns", "resultOutputColumns", ...
     "output_columns", "outputColumns", ...
     "table_columns", "tableColumns", ...
     "columns"]);
resultOutputColumnValues = stringListFieldValues( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["result_output_columns", "resultOutputColumns", ...
     "output_columns", "outputColumns", ...
     "table_columns", "tableColumns", ...
     "columns"]);
resultOutputUnits = firstUnitStructField( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["result_output_units", "resultOutputUnits", ...
     "output_units", "outputUnits", ...
     "table_units", "tableUnits", ...
     "column_units", "columnUnits", ...
     "units"]);
resultOutputUnitValues = unitStructFieldValues( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["result_output_units", "resultOutputUnits", ...
     "output_units", "outputUnits", ...
     "table_units", "tableUnits", ...
     "column_units", "columnUnits", ...
     "units"]);
expectedResultOutputUnits = normalizeUnitStruct(options.ExpectedResultOutputUnits);
resultOutputSchemaRequired = options.RequireResultOutputSchema || ...
    strlength(options.ExpectedResultOutputSchemaId) > 0 || ...
    ~isempty(options.ExpectedResultOutputColumns) || ...
    hasUnitFields(expectedResultOutputUnits);
physicsConventionSchemaId = firstStringField( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["physics_convention_schema_id", "physicsConventionSchemaId", ...
     "coupling_convention_schema_id", "couplingConventionSchemaId", ...
     "fem_bem_coupling_convention_schema_id", "femBemCouplingConventionSchemaId", ...
     "result_convention_schema_id", "resultConventionSchemaId"]);
physicsConventionSchemaIdValues = uniqueStringFields( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["physics_convention_schema_id", "physicsConventionSchemaId", ...
     "coupling_convention_schema_id", "couplingConventionSchemaId", ...
     "fem_bem_coupling_convention_schema_id", "femBemCouplingConventionSchemaId", ...
     "result_convention_schema_id", "resultConventionSchemaId"]);
physicsConventionSchemaRequired = options.RequirePhysicsConventionSchema || ...
    strlength(options.ExpectedPhysicsConventionSchemaId) > 0;
postprocessRowConventionSchemaId = firstStringField( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["postprocess_row_convention_schema_id", "postprocessRowConventionSchemaId", ...
     "fem_bem_postprocess_row_convention_schema_id", "femBemPostprocessRowConventionSchemaId", ...
     "trace_postprocess_row_convention_schema_id", "tracePostprocessRowConventionSchemaId"]);
postprocessRowConventionSchemaIdValues = uniqueStringFields( ...
    {artifact, execution, resultBlock, outputBlock}, ...
    ["postprocess_row_convention_schema_id", "postprocessRowConventionSchemaId", ...
     "fem_bem_postprocess_row_convention_schema_id", "femBemPostprocessRowConventionSchemaId", ...
     "trace_postprocess_row_convention_schema_id", "tracePostprocessRowConventionSchemaId"]);
postprocessRowConventionSchemaRequired = options.RequirePostprocessRowConventionSchema || ...
    strlength(options.ExpectedPostprocessRowConventionSchemaId) > 0;

checks = struct();
checks.schema_recorded = hasNonemptyField(artifact, "schema");
checks.created_at_utc_recorded = strlength(createdAt) > 0;
checks.created_at_utc_parseable = looksLikeIsoDateTime(createdAt);
checks.run_date_utc_recorded_when_required = ~options.RequireRunDate || runDateRecorded;
checks.run_date_utc_parseable = looksLikeIsoDateTime(runDate);
checks.expected_created_at_utc_matches = strlength(options.ExpectedCreatedAtUtc) == 0 || createdAt == options.ExpectedCreatedAtUtc;
checks.expected_run_date_utc_matches = strlength(options.ExpectedRunDateUtc) == 0 || runDate == options.ExpectedRunDateUtc;
checks.execution_session_id_consistent_when_present = numel(executionSessionIdValues) <= 1;
checks.execution_session_id_recorded_when_required = ...
    ~options.RequireExecutionSessionId || strlength(executionSessionId) > 0;
checks.execution_session_id_matches_expected = ...
    strlength(options.ExpectedExecutionSessionId) == 0 || ...
    (~isempty(executionSessionIdValues) && ...
    all(executionSessionIdValues == options.ExpectedExecutionSessionId));
checks.created_run_timestamp_skew_within_limit = isinf(options.MaxCreatedRunSkewSeconds) || ...
    (~isnan(createdRunSkewSeconds) && createdRunSkewSeconds <= options.MaxCreatedRunSkewSeconds);
checks.versions_mapping_recorded = isstruct(versions) && ~isempty(versionKeys);
checks.required_versions_recorded = isempty(missingVersions);
checks.radia_version_recorded = ~options.RequireRadiaVersion || radiaVersionRecorded;
checks.pass_recorded_true = artifactPasses(artifact);
checks.timing_breakdown_recorded = height(timingValues) >= options.MinTimingStages;
checks.timing_stage_count_reasonable = height(timingValues) <= options.MaxTimingStages;
checks.timing_values_nonnegative = timingValuesOk;
checks.required_timing_stages_recorded = isempty(missingTimingStages);
checks.dominant_timing_stages_available = dominantCount > 0;
checks.parameter_set_artifact_id_consistent_when_present = numel(parameterSetArtifactIdValues) <= 1;
checks.parameter_set_digest_consistent_when_present = numel(parameterSetDigestValues) <= 1;
checks.parameter_set_path_consistent_when_present = numel(parameterSetPathValues) <= 1;
checks.objective_observable_id_consistent_when_present = numel(objectiveObservableIdValues) <= 1;
checks.objective_observable_family_consistent_when_present = numel(objectiveObservableFamilyValues) <= 1;
checks.parameter_set_artifact_id_recorded_when_required = ...
    ~options.RequireParameterSetArtifact || strlength(parameterSetArtifactId) > 0;
checks.parameter_set_digest_recorded_when_required = ...
    ~options.RequireParameterSetArtifact || strlength(parameterSetDigest) > 0;
checks.parameter_set_path_recorded_when_required = ...
    ~options.RequireParameterSetArtifact || strlength(parameterSetPath) > 0;
checks.parameter_set_artifact_id_matches_expected = ...
    strlength(options.ExpectedParameterSetArtifactId) == 0 || ...
    (~isempty(parameterSetArtifactIdValues) && ...
    all(parameterSetArtifactIdValues == options.ExpectedParameterSetArtifactId));
checks.parameter_set_digest_matches_expected = ...
    strlength(options.ExpectedParameterSetDigest) == 0 || ...
    (~isempty(parameterSetDigestValues) && ...
    all(parameterSetDigestValues == options.ExpectedParameterSetDigest));
checks.parameter_set_path_matches_expected = ...
    strlength(options.ExpectedParameterSetPath) == 0 || ...
    (~isempty(parameterSetPathValues) && ...
    all(parameterSetPathValues == options.ExpectedParameterSetPath));
checks.objective_observable_id_matches_expected = ...
    strlength(options.ExpectedObjectiveObservableId) == 0 || ...
    (~isempty(objectiveObservableIdValues) && ...
    all(objectiveObservableIdValues == options.ExpectedObjectiveObservableId));
checks.objective_observable_family_matches_expected = ...
    strlength(options.ExpectedObjectiveObservableFamily) == 0 || ...
    (~isempty(objectiveObservableFamilyValues) && ...
    all(objectiveObservableFamilyValues == options.ExpectedObjectiveObservableFamily));
checks.result_output_schema_id_consistent_when_present = numel(resultOutputSchemaIdValues) <= 1;
checks.result_output_columns_consistent_when_present = allStringListValuesSame(resultOutputColumnValues);
checks.result_output_units_consistent_when_present = allUnitStructValuesSame(resultOutputUnitValues);
checks.result_output_schema_id_recorded_when_required = ...
    ~resultOutputSchemaRequired || strlength(resultOutputSchemaId) > 0;
checks.result_output_columns_recorded_when_required = ...
    ~resultOutputSchemaRequired || ~isempty(resultOutputColumns);
checks.result_output_units_recorded_when_required = ...
    ~resultOutputSchemaRequired || hasUnitFields(resultOutputUnits);
checks.result_output_schema_id_matches_expected = ...
    strlength(options.ExpectedResultOutputSchemaId) == 0 || ...
    resultOutputSchemaId == options.ExpectedResultOutputSchemaId;
checks.result_output_columns_match_expected = ...
    isempty(options.ExpectedResultOutputColumns) || ...
    isequal(resultOutputColumns(:), options.ExpectedResultOutputColumns(:));
checks.result_output_units_match_expected = ...
    ~hasUnitFields(expectedResultOutputUnits) || ...
    unitStructsEqual(resultOutputUnits, expectedResultOutputUnits);
checks.physics_convention_schema_id_consistent_when_present = ...
    numel(physicsConventionSchemaIdValues) <= 1;
checks.physics_convention_schema_id_recorded_when_required = ...
    ~physicsConventionSchemaRequired || strlength(physicsConventionSchemaId) > 0;
checks.physics_convention_schema_id_matches_expected = ...
    strlength(options.ExpectedPhysicsConventionSchemaId) == 0 || ...
    physicsConventionSchemaId == options.ExpectedPhysicsConventionSchemaId;
checks.postprocess_row_convention_schema_id_consistent_when_present = ...
    numel(postprocessRowConventionSchemaIdValues) <= 1;
checks.postprocess_row_convention_schema_id_recorded_when_required = ...
    ~postprocessRowConventionSchemaRequired || strlength(postprocessRowConventionSchemaId) > 0;
checks.postprocess_row_convention_schema_id_matches_expected = ...
    strlength(options.ExpectedPostprocessRowConventionSchemaId) == 0 || ...
    postprocessRowConventionSchemaId == options.ExpectedPostprocessRowConventionSchemaId;

status = "ok";
checkNames = fieldnames(checks);
for k = 1:numel(checkNames)
    if ~checks.(checkNames{k})
        status = "needs_attention";
        break
    end
end

report = struct();
report.policy = "matlab_result_manifest_gate";
report.status = status;
report.created_at_utc = createdAt;
report.run_date_utc = runDate;
report.expected_created_at_utc = options.ExpectedCreatedAtUtc;
report.expected_run_date_utc = options.ExpectedRunDateUtc;
report.require_run_date = options.RequireRunDate;
report.execution_session_id = executionSessionId;
report.execution_session_id_values = executionSessionIdValues;
report.expected_execution_session_id = options.ExpectedExecutionSessionId;
report.require_execution_session_id = options.RequireExecutionSessionId;
report.created_run_skew_s = createdRunSkewSeconds;
report.max_created_run_skew_s = options.MaxCreatedRunSkewSeconds;
report.required_versions = options.RequiredVersions;
report.missing_versions = missingVersions;
report.version_keys = versionKeys;
report.radia_version_recorded = radiaVersionRecorded;
report.timing_stage_count = height(timingValues);
report.max_timing_stages = options.MaxTimingStages;
report.required_timing_stages = options.RequiredTimingStages;
report.missing_timing_stages = missingTimingStages;
report.parameter_set_artifact_id = parameterSetArtifactId;
report.parameter_set_digest = parameterSetDigest;
report.parameter_set_path = parameterSetPath;
report.objective_observable_id = objectiveObservableId;
report.objective_observable_family = objectiveObservableFamily;
report.parameter_set_artifact_id_values = parameterSetArtifactIdValues;
report.parameter_set_digest_values = parameterSetDigestValues;
report.parameter_set_path_values = parameterSetPathValues;
report.objective_observable_id_values = objectiveObservableIdValues;
report.objective_observable_family_values = objectiveObservableFamilyValues;
report.require_parameter_set_artifact = options.RequireParameterSetArtifact;
report.result_output_schema_id = resultOutputSchemaId;
report.result_output_schema_id_values = resultOutputSchemaIdValues;
report.result_output_columns = resultOutputColumns;
report.result_output_units = resultOutputUnits;
report.expected_result_output_schema_id = options.ExpectedResultOutputSchemaId;
report.expected_result_output_columns = options.ExpectedResultOutputColumns;
report.expected_result_output_units = expectedResultOutputUnits;
report.require_result_output_schema = options.RequireResultOutputSchema;
report.physics_convention_schema_id = physicsConventionSchemaId;
report.physics_convention_schema_id_values = physicsConventionSchemaIdValues;
report.expected_physics_convention_schema_id = options.ExpectedPhysicsConventionSchemaId;
report.require_physics_convention_schema = options.RequirePhysicsConventionSchema;
report.postprocess_row_convention_schema_id = postprocessRowConventionSchemaId;
report.postprocess_row_convention_schema_id_values = postprocessRowConventionSchemaIdValues;
report.expected_postprocess_row_convention_schema_id = options.ExpectedPostprocessRowConventionSchemaId;
report.require_postprocess_row_convention_schema = options.RequirePostprocessRowConventionSchema;
report.timing_breakdown_s = timingValues;
report.dominant_timing_stages = dominantTiming;
report.total_recorded_timing_s = sum(timingValues.seconds);
report.checks = checks;
report.notes = [
    "Record run date and versions before importing JSON into notebooks."
    "Keep created_at_utc and execution.run_date_utc close enough to describe the same executed result."
    "Record the execution session id when a notebook result depends on an already-running MATLAB session."
    "Keep timing to the heaviest stages, normally setup/assembly/solve/postprocess."
    "Record parameter-set identity and objective observable identity before a notebook reuses JSON values as defaults or optimization inputs."
    "Record result output schema id, columns, and units before a notebook imports saved JSON/table rows."
    "Record a physics or FEM/BEM coupling convention schema id before reusing values whose meaning depends on normals, trace rows, spaces, sources, or kernels."
    "Record a postprocess row convention schema id before reusing trace least-squares rows, residual curves, scalar objectives, or optimizer-facing table reductions."
    "Store compact JSON sidecars, then render them into executed ipynb panels or docs."
];
end

function artifact = loadArtifact(result)
if ischar(result) || (isstring(result) && isscalar(result))
    text = fileread(char(result));
    artifact = jsondecode(text);
elseif isstruct(result)
    artifact = result;
else
    error("result_manifest_gate:InvalidInput", "result must be a struct or JSON file path")
end
if ~isscalar(artifact)
    error("result_manifest_gate:InvalidInput", "result artifact must be a scalar struct")
end
end

function value = getStructField(s, fieldName, defaultValue)
name = char(fieldName);
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = defaultValue;
end
end

function tf = hasNonemptyField(s, fieldName)
name = matlab.lang.makeValidName(char(fieldName));
tf = isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && strlength(string(s.(name))) > 0;
end

function tf = hasAnyNonemptyField(s, fieldNames)
tf = false;
for key = string(fieldNames)
    if hasNonemptyField(s, key)
        tf = true;
        return
    end
end
end

function value = firstStringField(records, fieldNames)
value = "";
for r = 1:numel(records)
    record = records{r};
    for key = string(fieldNames)
        if hasNonemptyField(record, key)
            value = string(record.(matlab.lang.makeValidName(char(key))));
            return
        end
    end
end
end

function values = uniqueStringFields(records, fieldNames)
values = strings(0, 1);
for r = 1:numel(records)
    record = records{r};
    for key = string(fieldNames)
        if hasNonemptyField(record, key)
            values(end + 1, 1) = string(record.(matlab.lang.makeValidName(char(key)))); %#ok<AGROW>
        end
    end
end
values = unique(values);
end

function value = firstStringListField(records, fieldNames)
value = strings(1, 0);
for r = 1:numel(records)
    record = records{r};
    for key = string(fieldNames)
        name = matlab.lang.makeValidName(char(key));
        if isstruct(record) && isfield(record, name) && ~isempty(record.(name))
            value = asStringList(record.(name));
            if ~isempty(value)
                return
            end
        end
    end
end
end

function values = stringListFieldValues(records, fieldNames)
values = {};
for r = 1:numel(records)
    record = records{r};
    for key = string(fieldNames)
        name = matlab.lang.makeValidName(char(key));
        if isstruct(record) && isfield(record, name) && ~isempty(record.(name))
            item = asStringList(record.(name));
            if ~isempty(item)
                values{end + 1} = item; %#ok<AGROW>
            end
        end
    end
end
end

function value = asStringList(raw)
if isstring(raw)
    value = raw(:).';
elseif iscell(raw)
    value = string(raw(:)).';
elseif ischar(raw)
    text = string(raw);
    if contains(text, ",") || contains(text, ";")
        value = strip(split(replace(text, ";", ","), ",")).';
    else
        value = text;
    end
elseif isnumeric(raw) || islogical(raw)
    value = string(raw(:)).';
else
    value = string(raw);
end
value = value(strlength(value) > 0);
end

function tf = allStringListValuesSame(values)
tf = true;
if numel(values) <= 1
    return
end
first = values{1};
for k = 2:numel(values)
    if ~isequal(first(:), values{k}(:))
        tf = false;
        return
    end
end
end

function value = firstUnitStructField(records, fieldNames)
value = struct();
for r = 1:numel(records)
    record = records{r};
    for key = string(fieldNames)
        name = matlab.lang.makeValidName(char(key));
        if isstruct(record) && isfield(record, name) && ~isempty(record.(name))
            item = normalizeUnitStruct(record.(name));
            if hasUnitFields(item)
                value = item;
                return
            end
        end
    end
end
end

function values = unitStructFieldValues(records, fieldNames)
values = {};
for r = 1:numel(records)
    record = records{r};
    for key = string(fieldNames)
        name = matlab.lang.makeValidName(char(key));
        if isstruct(record) && isfield(record, name) && ~isempty(record.(name))
            item = normalizeUnitStruct(record.(name));
            if hasUnitFields(item)
                values{end + 1} = item; %#ok<AGROW>
            end
        end
    end
end
end

function value = normalizeUnitStruct(raw)
value = struct();
if ~isstruct(raw)
    return
end
names = fieldnames(raw);
for k = 1:numel(names)
    name = matlab.lang.makeValidName(names{k});
    unit = raw.(names{k});
    if ~isempty(name)
        value.(name) = string(unit);
    end
end
end

function tf = hasUnitFields(value)
tf = isstruct(value) && ~isempty(fieldnames(value));
end

function tf = unitStructsEqual(a, b)
a = normalizeUnitStruct(a);
b = normalizeUnitStruct(b);
fieldsA = sort(string(fieldnames(a)));
fieldsB = sort(string(fieldnames(b)));
if ~isequal(fieldsA, fieldsB)
    tf = false;
    return
end
tf = true;
for field = fieldsA.'
    name = char(field);
    if string(a.(name)) ~= string(b.(name))
        tf = false;
        return
    end
end
end

function tf = allUnitStructValuesSame(values)
tf = true;
if numel(values) <= 1
    return
end
first = values{1};
for k = 2:numel(values)
    if ~unitStructsEqual(first, values{k})
        tf = false;
        return
    end
end
end

function tf = looksLikeIsoDateTime(value)
text = char(string(value));
tf = ~isempty(regexp(text, "^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", "once"));
end

function seconds = isoDateTimeSeconds(value)
text = char(string(value));
tokens = regexp(text, "^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})", "tokens", "once");
if isempty(tokens)
    seconds = NaN;
    return
end
parts = str2double(tokens);
stamp = datetime(parts(1), parts(2), parts(3), parts(4), parts(5), parts(6), "TimeZone", "UTC");
seconds = posixtime(stamp);
end

function tf = artifactPasses(artifact)
if ~isfield(artifact, "pass")
    tf = false;
    return
end
value = artifact.pass;
if islogical(value)
    tf = isscalar(value) && value;
elseif isnumeric(value)
    tf = isscalar(value) && value ~= 0;
elseif ischar(value) || isstring(value)
    tf = any(lower(string(value)) == ["true", "pass", "passed", "ok"]);
else
    tf = false;
end
end

function [timingValues, valuesOk] = timingStructToTable(timing)
valuesOk = true;
stages = strings(0, 1);
seconds = zeros(0, 1);
if ~isstruct(timing)
    valuesOk = false;
    timingValues = table(stages, seconds, 'VariableNames', {'stage', 'seconds'});
    return
end
fields = fieldnames(timing);
for k = 1:numel(fields)
    value = timing.(fields{k});
    if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
        valuesOk = false;
        continue
    end
    stages(end + 1, 1) = string(fields{k}); %#ok<AGROW>
    seconds(end + 1, 1) = double(value); %#ok<AGROW>
end
timingValues = table(stages, seconds, 'VariableNames', {'stage', 'seconds'});
end
