function check_result_manifest_file(manifestPath, requireRadiaVersion, requireResultOutputSchema, expectedResultOutputSchemaId)
%CHECK_RESULT_MANIFEST_FILE Print a compact JSON report for a result manifest.
%
% This MCP-facing entry point is intended for MathWorks MATLAB MCP Server custom-tool
% extension files.  The official extension format accepts scalar arguments, so
% larger result data should be passed by file path.

arguments
    manifestPath (1,1) string
    requireRadiaVersion (1,1) logical = true
    requireResultOutputSchema (1,1) logical = false
    expectedResultOutputSchemaId (1,1) string = ""
end

report = acoustic_fembem.result_manifest_gate( ...
    manifestPath, ...
    RequireRadiaVersion=requireRadiaVersion, ...
    RequireResultOutputSchema=requireResultOutputSchema, ...
    ExpectedResultOutputSchemaId=expectedResultOutputSchemaId);

summary = struct();
summary.tool = "acoustic_fembem_check_result_manifest_file";
summary.status = report.status;
summary.ok = report.status == "ok";
summary.policy = report.policy;
summary.manifest_path = manifestPath;
summary.created_at_utc = report.created_at_utc;
summary.run_date_utc = report.run_date_utc;
summary.execution_session_id = report.execution_session_id;
summary.result_output_schema_id = report.result_output_schema_id;
summary.physics_convention_schema_id = report.physics_convention_schema_id;
summary.postprocess_row_convention_schema_id = report.postprocess_row_convention_schema_id;
summary.dominant_timing_stages = tableToStructArray(report.dominant_timing_stages);
summary.failed_checks = failedCheckNames(report.checks);

disp(jsonencode(summary));

if summary.ok
    return
end

error("acoustic_fembem:ResultManifestNeedsAttention", ...
    "Result manifest needs attention. Failed checks: %s", strjoin(summary.failed_checks, ", "));
end

function names = failedCheckNames(checks)
names = strings(1, 0);
checkNames = string(fieldnames(checks));
for k = 1:numel(checkNames)
    checkName = checkNames(k);
    if ~checks.(checkName)
        names(end + 1) = checkName; %#ok<AGROW>
    end
end
end

function rows = tableToStructArray(value)
if istable(value)
    rows = table2struct(value);
else
    rows = struct([]);
end
end
