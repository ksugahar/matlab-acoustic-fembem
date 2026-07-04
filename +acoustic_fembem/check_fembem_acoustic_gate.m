function check_fembem_acoustic_gate(kind, wavenumber, quadratureOrder, tolerance)
%CHECK_FEMBEM_ACOUSTIC_GATE Print a compact JSON verdict for an acoustic gate.
%
% MCP-facing entry point for the MathWorks MATLAB MCP Server custom-tool extension:
% runs the integrated Gypsilab acoustic solver against the analytic
% partial-wave series and prints a scalar-friendly JSON verdict. The
% official extension format accepts scalar arguments, so knobs are scalars;
% pass tolerance = NaN (or a nonpositive value) to use the per-kind default.

arguments
    kind (1,1) string
    wavenumber (1,1) double = 2.0
    quadratureOrder (1,1) double = 7
    tolerance (1,1) double = -1
end

args = {"Wavenumber", wavenumber, "QuadratureOrder", quadratureOrder};
if tolerance > 0
    args = [args, {"Tolerance", tolerance}];
end
report = acoustic_fembem.fembem_acoustic_gate(kind, args{:});

summary = struct();
summary.tool = report.tool;
summary.status = report.status;
summary.ok = report.pass;
summary.kind = report.kind;
summary.fixture = report.fixture;
summary.wavenumber = report.wavenumber;
summary.quadrature_order = report.quadrature_order;
summary.reference = report.reference;
summary.relative_error = report.relative_error;
summary.tolerance = report.tolerance;

disp(jsonencode(summary));

if summary.ok
    return
end
error("acoustic_fembem:FembemAcousticGateNeedsAttention", ...
    "Acoustic gate '%s' rel err %.3e exceeds tolerance %.3e.", ...
    report.kind, report.relative_error, report.tolerance);
end
