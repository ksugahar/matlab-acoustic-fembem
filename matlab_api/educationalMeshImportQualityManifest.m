function report = educationalMeshImportQualityManifest(surfaceElementTypes, volumeElementTypes, varargin)
%EDUCATIONALMESHIMPORTQUALITYMANIFEST Readable mesh handoff preflight.
%
% This helper is a teaching manifest before a third-party mesh is promoted to
% the MATLAB Gypsilab/Lukas-readable lane.  It deliberately keeps the MATLAB
% parser policy strict: first-order boundary triangles and volume tetrahedra
% only.  High-order, quad, hex, prism, wedge, or pyramid evidence can be logged
% here, but it is not silently converted into the .vol tri/tet intake.

if nargin < 1 || isempty(surfaceElementTypes)
    surfaceElementTypes = "tri";
end
if nargin < 2 || isempty(volumeElementTypes)
    volumeElementTypes = "tet";
end

p = inputParser;
p.FunctionName = "educationalMeshImportQualityManifest";
addParameter(p, "Source", "third_party_mesh", @(x) isstring(x) || ischar(x));
addParameter(p, "Format", ".vol", @(x) isstring(x) || ischar(x));
addParameter(p, "Order", 1, @(x) isscalar(x) && isfinite(x));
addParameter(p, "MinScaledJacobianBefore", NaN, @(x) isscalar(x) && (isnan(x) || isfinite(x)));
addParameter(p, "MinScaledJacobianAfter", 1.0, @(x) isscalar(x) && isfinite(x));
addParameter(p, "MinScaledJacobianThreshold", 0.1, @(x) isscalar(x) && isfinite(x));
addParameter(p, "NegativeJacobianCountBefore", 0, @(x) isscalar(x) && isfinite(x));
addParameter(p, "NegativeJacobianCountAfter", 0, @(x) isscalar(x) && isfinite(x));
addParameter(p, "CadConnectivityRecorded", true, @(x) islogical(x) && isscalar(x));
addParameter(p, "CadComplianceRecorded", true, @(x) islogical(x) && isscalar(x));
addParameter(p, "BoundaryConformityTolerance", 1e-8, @(x) isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, "MaxBoundaryDistance", 0.0, @(x) isscalar(x) && isfinite(x) && x >= 0);
parse(p, varargin{:});
opts = p.Results;

surface = normalizeElementTypes(surfaceElementTypes);
volume = normalizeElementTypes(volumeElementTypes);
order = double(opts.Order);
negativeBefore = double(opts.NegativeJacobianCountBefore);
negativeAfter = double(opts.NegativeJacobianCountAfter);
minJacBefore = double(opts.MinScaledJacobianBefore);
minJacAfter = double(opts.MinScaledJacobianAfter);
minJacThreshold = double(opts.MinScaledJacobianThreshold);

checks = struct();
checks.sourceRecorded = strlength(string(opts.Source)) > 0;
checks.formatRecorded = strlength(string(opts.Format)) > 0;
checks.surfaceTrianglesOnly = ~isempty(surface) && all(surface == "tri");
checks.volumeTetrahedraOnly = ~isempty(volume) && all(volume == "tet");
checks.firstOrderOnly = isfinite(order) && order == round(order) && order == 1;
checks.noNegativeJacobianAfter = negativeAfter == 0;
checks.negativeJacobianCountImproved = negativeAfter <= negativeBefore;
checks.minScaledJacobianRecorded = isfinite(minJacAfter);
checks.minScaledJacobianAboveThreshold = minJacAfter >= minJacThreshold;
checks.cadConnectivityRecorded = opts.CadConnectivityRecorded;
checks.cadComplianceRecorded = opts.CadComplianceRecorded;
checks.boundaryConformityWithinTolerance = opts.MaxBoundaryDistance <= opts.BoundaryConformityTolerance;
checks.matlabVolParserPolicyHonored = checks.surfaceTrianglesOnly ...
    && checks.volumeTetrahedraOnly ...
    && checks.firstOrderOnly;

issues = strings(0, 1);
names = string(fieldnames(checks));
for k = 1:numel(names)
    if ~checks.(names(k))
        issues(end + 1, 1) = names(k); %#ok<AGROW>
    end
end

if isfinite(minJacBefore)
    jacobianImprovement = minJacAfter - minJacBefore;
else
    jacobianImprovement = NaN;
end

report = struct();
report.kind = "educational_mesh_import_quality_manifest";
report.policy = "readable_mesh_import_manifest_first_order_tri_tet_only";
report.source = string(opts.Source);
report.format = string(opts.Format);
report.surfaceElementTypes = surface;
report.volumeElementTypes = volume;
report.order = order;
report.minScaledJacobianBefore = minJacBefore;
report.minScaledJacobianAfter = minJacAfter;
report.minScaledJacobianThreshold = minJacThreshold;
report.minScaledJacobianImprovement = jacobianImprovement;
report.negativeJacobianCountBefore = negativeBefore;
report.negativeJacobianCountAfter = negativeAfter;
report.boundaryConformityTolerance = opts.BoundaryConformityTolerance;
report.maxBoundaryDistance = opts.MaxBoundaryDistance;
report.checks = checks;
report.issues = issues;
report.pass = isempty(issues);
if report.pass
    report.status = "ok";
else
    report.status = "needs_attention";
end
end


function out = normalizeElementTypes(values)
out = lower(strtrim(string(values(:))));
out(out == "triangle") = "tri";
out(out == "tri3") = "tri";
out(out == "tetra") = "tet";
out(out == "tetrahedron") = "tet";
out(out == "tet4") = "tet";
end
