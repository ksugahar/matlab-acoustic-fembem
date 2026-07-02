function report = femBemTraceMassReport(input, femValues, tol)
%FEMBEMTRACEMASSREPORT Check P1 FEM trace and BEM surface mass.
%
% The readable identity is
%
%   (T u)' M_b (T u) == u' T' M_b T u,
%
% where T maps volume H1/P1 nodes to compact boundary P1 nodes and M_b is the
% boundary surface mass matrix.  This is a small teaching gate for FEM/BEM
% coupling before adding singular BEM kernels.

if nargin < 3
    tol = 1e-12;
end

if ischar(input) || isstring(input)
    model = FemBemModel(string(input));
elseif isa(input, "FemBemModel")
    model = input;
else
    error("femBemTraceMassReport:input", ...
        "Input must be a .vol path or a FemBemModel.");
end

if isempty(model.operators)
    model = model.assemble();
end

T = model.operators.trace.matrix;
Mb = model.operators.bem.surfaceMass;
nFem = size(T, 2);
nBem = size(T, 1);
meshId = model.mesh.meshId;
surfaceMeshId = model.operators.trace.surfaceMeshId;
traceArtifactId = model.operators.trace.artifactId;

if nargin < 2 || isempty(femValues)
    femValues = (1:nFem).';
else
    femValues = femValues(:);
end

if numel(femValues) ~= nFem
    error("femBemTraceMassReport:size", ...
        "femValues must have one entry per FEM node.");
end

boundaryValues = T * femValues;
massBoundaryValues = Mb * boundaryValues;
embeddedBoundaryAction = T' * massBoundaryValues;
boundaryEnergy = boundaryValues' * massBoundaryValues;
embeddedEnergy = femValues' * embeddedBoundaryAction;
onesFem = ones(nFem, 1);
onesBem = ones(nBem, 1);
tracedConstant = T * onesFem;
surfaceAreaFromMass = onesBem' * Mb * onesBem;
surfaceMassSymmetryError = full(max(max(abs(Mb - Mb.'))));
constantTraceError = max(abs(tracedConstant - 1));
interiorNodeIds = find(full(sum(abs(T), 1)) == 0).';
if isempty(interiorNodeIds)
    interiorBoundaryAction = 0;
else
    interiorBoundaryAction = max(abs(embeddedBoundaryAction(interiorNodeIds)));
end
energyAbsError = abs(boundaryEnergy - embeddedEnergy);
energyRelError = energyAbsError / max([abs(boundaryEnergy), abs(embeddedEnergy), realmin]);

checks = struct();
checks.surfaceMassSymmetric = surfaceMassSymmetryError <= tol;
checks.constantTraceOk = constantTraceError <= tol;
checks.energyIdentityOk = energyAbsError <= tol || energyRelError <= tol;
checks.interiorNodesHaveZeroBoundaryAction = interiorBoundaryAction <= tol;
checks.surfaceAreaPositive = surfaceAreaFromMass > 0;
checks.meshIdRecorded = meshId ~= "";
checks.surfaceMeshIdRecorded = surfaceMeshId ~= "";
checks.traceArtifactIdRecorded = traceArtifactId ~= "";

report = struct();
report.kind = "fem_bem_trace_mass_report";
report.policy = "readable_first_order_tri_tet_trace_surface_mass_gate";
report.meshId = meshId;
report.surfaceMeshId = surfaceMeshId;
report.traceArtifactId = traceArtifactId;
report.traceShape = [nBem, nFem];
report.surfaceMassShape = size(Mb);
report.interiorNodeIds = interiorNodeIds;
report.boundaryValues = boundaryValues;
report.embeddedBoundaryAction = embeddedBoundaryAction;
report.boundaryEnergy = boundaryEnergy;
report.embeddedEnergy = embeddedEnergy;
report.energyAbsError = energyAbsError;
report.energyRelError = energyRelError;
report.surfaceAreaFromMass = surfaceAreaFromMass;
report.surfaceMassSymmetryError = surfaceMassSymmetryError;
report.constantTraceError = constantTraceError;
report.interiorBoundaryAction = interiorBoundaryAction;
report.checks = checks;
report.status = string(statusFromChecks(checks));
end


function status = statusFromChecks(checks)
values = structfun(@(value) logical(value), checks);
if all(values)
    status = "ok";
else
    status = "needs_attention";
end
end
