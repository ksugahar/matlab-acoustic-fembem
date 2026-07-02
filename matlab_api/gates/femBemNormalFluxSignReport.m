function report = femBemNormalFluxSignReport(model, gradientVector, options)
%FEMBEMNORMALFLUXSIGNREPORT Check .vol outward-normal flux signs.
%
% The .vol reader may store boundary triangles inward or outward. This helper
% keeps the sign correction visible before a FEM/BEM normal derivative, flux
% integral, or surface source row is used in a notebook.

arguments
    model (1,1) FemBemModel
    gradientVector (1,3) double = [1 2 3]
    options.AbsTol (1,1) double = 1e-12
    options.NormalConvention (1,1) string = "outward_from_volume"
    options.ExpectedNormalConvention (1,1) string = "outward_from_volume"
    options.ResultArtifactId (1,1) string = "matlab_fem_bem_normal_flux_sign_report_v1"
    options.ExpectedResultArtifactId (1,1) string = ""
    options.ResultDigest (1,1) string = ""
    options.ExpectedResultDigest (1,1) string = ""
end

orientation = model.surface.orientation;
rows = orientation.rows;
signs = double(orientation.triangleOrientationSignsToOutward(:));
nRows = numel(rows);
storedFlux = zeros(nRows, 1);
outwardReferenceFlux = zeros(nRows, 1);

for k = 1:nRows
    storedAreaVector = double(rows(k).storedAreaVector);
    outwardAreaVector = double(rows(k).outwardAreaVector);
    storedFlux(k) = dot(gradientVector, storedAreaVector);
    outwardReferenceFlux(k) = dot(gradientVector, outwardAreaVector);
end

orientationCorrectedFlux = signs .* storedFlux;
localAbsError = abs(orientationCorrectedFlux - outwardReferenceFlux);
closedSurfaceFluxSum = sum(orientationCorrectedFlux);
expectedClosedSurfaceFluxSum = 0.0;
tolerance = options.AbsTol;
normalConvention = lower(strrep(strtrim(options.NormalConvention), "-", "_"));
expectedNormalConvention = lower(strrep(strtrim(options.ExpectedNormalConvention), "-", "_"));
resultArtifactId = strtrim(options.ResultArtifactId);
expectedResultArtifactId = defaultExpected(options.ExpectedResultArtifactId, resultArtifactId);
resultDigest = strtrim(options.ResultDigest);
if strlength(resultDigest) == 0
    resultDigest = localNormalFluxDigest( ...
        signs, ...
        storedFlux, ...
        orientationCorrectedFlux, ...
        outwardReferenceFlux, ...
        normalConvention);
end
expectedResultDigest = defaultExpected(options.ExpectedResultDigest, resultDigest);

checks = struct();
checks.orientationSignsRecorded = numel(signs) == nRows && nRows > 0;
checks.orientationSignsAreUnit = ~isempty(signs) && all(abs(signs) == 1);
checks.normalConventionRecorded = strlength(normalConvention) > 0;
checks.normalConventionMatchesExpected = normalConvention == expectedNormalConvention;
checks.resultArtifactIdRecorded = strlength(resultArtifactId) > 0;
checks.resultArtifactIdMatchesExpected = resultArtifactId == expectedResultArtifactId;
checks.resultDigestRecorded = strlength(resultDigest) > 0;
checks.resultDigestMatchesExpected = resultDigest == expectedResultDigest;
checks.correctedFluxMatchesOutwardReference = max(localAbsError, [], "omitnan") <= tolerance;
checks.closedSurfaceFluxBalanceOk = abs(closedSurfaceFluxSum - expectedClosedSurfaceFluxSum) <= tolerance;

if all(struct2array(checks))
    status = "ok";
else
    status = "needs_attention";
end

report = struct();
report.policy = "readable_fem_bem_normal_flux_orientation_gate";
report.status = status;
report.resultArtifactId = resultArtifactId;
report.expectedResultArtifactId = expectedResultArtifactId;
report.resultDigest = resultDigest;
report.expectedResultDigest = expectedResultDigest;
report.gradientVector = gradientVector;
report.normalConvention = normalConvention;
report.expectedNormalConvention = expectedNormalConvention;
report.triangleOrientationSignsToOutward = signs;
report.storedNormalFlux = storedFlux;
report.orientationCorrectedNormalFlux = orientationCorrectedFlux;
report.outwardNormalFluxReference = outwardReferenceFlux;
report.localAbsError = localAbsError;
report.maxAbsError = max(localAbsError, [], "omitnan");
report.closedSurfaceFluxSum = closedSurfaceFluxSum;
report.expectedClosedSurfaceFluxSum = expectedClosedSurfaceFluxSum;
report.closedSurfaceFluxResidual = abs(closedSurfaceFluxSum - expectedClosedSurfaceFluxSum);
report.tolerance = tolerance;
report.checks = checks;
report.version = version;
report.runDate = string(datetime("now", "TimeZone", "local", "Format", "yyyy-MM-dd'T'HH:mm:ssXXX"));
end


function expected = defaultExpected(value, fallback)
if strlength(value) == 0
    expected = fallback;
else
    expected = value;
end
end


function digest = localNormalFluxDigest(signs, storedFlux, correctedFlux, referenceFlux, normalConvention)
lines = [
    "normal_flux_sign_report_v1"
    "normal_convention=" + string(normalConvention)
    ];
for k = 1:numel(signs)
    lines(end + 1, 1) = "row=" + string(k) + ...
        ",sign=" + string(signs(k)) + ...
        ",stored=" + sprintf("%.17g", storedFlux(k)) + ...
        ",corrected=" + sprintf("%.17g", correctedFlux(k)) + ...
        ",reference=" + sprintf("%.17g", referenceFlux(k)); %#ok<AGROW>
end
payload = join(lines, newline);
md = javaMethod("getInstance", "java.security.MessageDigest", "SHA-256");
md.update(uint8(char(payload)));
hash = typecast(md.digest(), "uint8");
hex = lower(reshape(dec2hex(hash, 2).', 1, []));
digest = "sha256:" + string(hex);
end
