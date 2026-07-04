function check_fembem_crossval_gate(kind, fixture, tolerance, writeLog)
%CHECK_FEMBEM_CROSSVAL_GATE Print a JSON verdict for a radia-ngsolve gate.
%
% MCP-facing entry point for cross-validation between the readable MATLAB
% FEM/BEM repository and the radia-ngsolve/NGSolve validation ladder.

arguments
    kind (1,1) string
    fixture (1,1) string = ""
    tolerance (1,1) double = -1
    writeLog (1,1) logical = false
end

args = {"Fixture", fixture, "WriteLog", writeLog};
if tolerance > 0
    args = [args, {"Tolerance", tolerance}];
end
report = acoustic_fembem.fembem_crossval_gate(kind, args{:});

summary = struct();
summary.tool = report.tool;
summary.status = report.status;
summary.ok = report.pass;
summary.kind = report.kind;
summary.fixture = report.fixture;
summary.reference = report.reference;
summary.reference_family = report.reference_family;
summary.input_format = report.input_format;
summary.relative_error = report.relative_error;
summary.tolerance = report.tolerance;
summary.duration_s = report.duration_s;
if isfield(report, "ngsolve_version")
    summary.ngsolve_version = report.ngsolve_version;
end
if isfield(report, "num_cases")
    summary.num_cases = report.num_cases;
end
if isfield(report, "num_passed")
    summary.num_passed = report.num_passed;
end

disp(jsonencode(summary));

if summary.ok
    return
end
error("acoustic_fembem:FembemCrossvalGateNeedsAttention", ...
    "Cross-validation gate '%s' rel err %.3e exceeds tolerance %.3e.", ...
    report.kind, report.relative_error, report.tolerance);
end
