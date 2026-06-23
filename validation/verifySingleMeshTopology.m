function result = verifySingleMeshTopology(item)
%VERIFYSINGLEMESHTOPOLOGY Compare MATLAB .vol topology with NGSolve.

arguments
    item (1,1) struct
end

matlabSummary = matlabMeshTopologySummary(item.volFile);
ngsolveSummary = ngsolveVolSummary(item.volFile);
failures = strings(0, 1);

if matlabSummary.ok ~= item.expectOk
    failures(end + 1, 1) = "MATLAB ok flag mismatch";
end
if item.expectOk && matlabSummary.ok && ngsolveSummary.ok
    failures = compareField(failures, matlabSummary, ngsolveSummary, "points");
    failures = compareField(failures, matlabSummary, ngsolveSummary, "triangles");
    failures = compareField(failures, matlabSummary, ngsolveSummary, "tets");
    failures = compareField(failures, matlabSummary, ngsolveSummary, "materials");
    failures = compareField(failures, matlabSummary, ngsolveSummary, "boundaries");
    failures = compareField(failures, matlabSummary, ngsolveSummary, "hcurlEdges");
    if ~isequal(sort(matlabSummary.materialNames), sort(string(ngsolveSummary.materialNames)))
        failures(end + 1, 1) = "material names mismatch";
    end
    if ~isequal(sort(matlabSummary.boundaryNames), sort(string(ngsolveSummary.boundaryNames)))
        failures(end + 1, 1) = "boundary names mismatch";
    end
elseif item.expectOk && ~ngsolveSummary.ok
    failures(end + 1, 1) = "NGSolve rejected an expected-ok tri/tet case";
elseif ~item.expectOk
    % NGSolve itself may accept broader Netgen element families. The education
    % solver policy is stricter: first-order tri/tet only.
    if matlabSummary.ok
        failures(end + 1, 1) = "MATLAB accepted a policy-rejected mesh";
    end
end

result = struct();
result.id = item.id;
result.title = item.title;
result.volFile = item.volFile;
result.expectOk = item.expectOk;
result.matlab = matlabSummary;
result.ngsolve = ngsolveSummary;
result.failures = failures;
result.passed = isempty(failures);
end


function failures = compareField(failures, matlabSummary, ngsolveSummary, field)
if matlabSummary.(field) ~= ngsolveSummary.(field)
    failures(end + 1, 1) = field + " mismatch";
end
end
