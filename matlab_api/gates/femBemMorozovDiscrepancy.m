function choice = femBemMorozovDiscrepancy(model, target, alphas, noiseNorm, options)
%femBemMorozovDiscrepancy Pick Tikhonov alpha from a noise norm.
%
% The path is computed by femBemTikhonovPath.  Morozov's readable
% teaching rule is then: choose the row whose weighted trace residual is
% closest to the estimated noise norm.  This pairs with the L-curve helper but
% keeps the noise-assumption choice explicit for students.

arguments
    model
    target double
    alphas double {mustBeNonnegative}
    noiseNorm (1,1) double {mustBePositive}
    options.Initial double = []
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.GradientTolerance (1,1) double {mustBePositive} = 1e-7
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-12
end

alphaValues = alphas(:);
if isempty(alphaValues)
    error("femBemMorozovDiscrepancy:alphas", ...
        "At least one Tikhonov weight is required.");
end
if any(diff(alphaValues) < 0)
    error("femBemMorozovDiscrepancy:alphas", ...
        "Tikhonov weights must be sorted in nondecreasing order.");
end

path = femBemTikhonovPath(model, target, alphaValues, ...
    "Initial", options.Initial, ...
    "FiniteDifferenceStep", options.FiniteDifferenceStep, ...
    "GradientTolerance", options.GradientTolerance);

residuals = path.weightedTraceResiduals(:);
errors = abs(residuals - noiseNorm);
[selectedAbsError, selectedIndex] = min(errors);

lowerIndex = find(residuals <= noiseNorm, 1, "last");
upperIndex = find(residuals >= noiseNorm, 1, "first");
noiseBracketed = ~isempty(lowerIndex) && ~isempty(upperIndex);

choice = struct();
choice.kind = "fem_bem_morozov_discrepancy";
choice.policy = "readable_morozov_discrepancy_from_trace_regularization_path";
choice.path = path;
choice.alphas = alphaValues;
choice.noiseNorm = noiseNorm;
choice.residualNorms = residuals;
choice.selectedIndex = selectedIndex;
choice.selectedAlpha = alphaValues(selectedIndex);
choice.selectedResidualNorm = residuals(selectedIndex);
choice.selectedAbsError = selectedAbsError;
choice.lowerBracketIndex = lowerIndex;
choice.upperBracketIndex = upperIndex;
choice.lowerBracketAlpha = bracketAlpha(alphaValues, lowerIndex);
choice.upperBracketAlpha = bracketAlpha(alphaValues, upperIndex);
choice.checks = struct( ...
    "pathOk", path.status == "ok", ...
    "alphasSorted", all(diff(alphaValues) >= 0), ...
    "residualsNondecreasing", all(diff(residuals) >= -options.Tolerance), ...
    "noiseBracketed", noiseBracketed, ...
    "selectedMinimizesDiscrepancy", selectedAbsError == min(errors));
choice.status = string(ifelse(all(struct2array(choice.checks)), "ok", "needs_attention"));
end


function value = bracketAlpha(alphas, index)
if isempty(index)
    value = NaN;
else
    value = alphas(index);
end
end


function value = ifelse(condition, a, b)
if condition
    value = a;
else
    value = b;
end
end
