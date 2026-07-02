function corner = femBemLcurveCorner(model, target, alphas, options)
%femBemLcurveCorner Pick a readable Tikhonov L-curve corner.
%
% The path is computed by femBemTikhonovPath.  This helper adds
% only the student-visible L-curve rule: take log(trace residual) and
% log(solution norm), compute a three-point discrete curvature, and select
% the largest interior curvature.  Endpoints are never selected.

arguments
    model
    target double
    alphas double {mustBePositive}
    options.Initial double = []
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.GradientTolerance (1,1) double {mustBePositive} = 1e-7
    options.LogFloor (1,1) double {mustBePositive} = 1e-300
end

alphaValues = alphas(:);
if numel(alphaValues) < 3
    error("femBemLcurveCorner:alphas", ...
        "At least three positive Tikhonov weights are required.");
end
if any(diff(alphaValues) <= 0)
    error("femBemLcurveCorner:alphas", ...
        "Tikhonov weights must be strictly increasing.");
end

path = femBemTikhonovPath(model, target, alphaValues, ...
    "Initial", options.Initial, ...
    "FiniteDifferenceStep", options.FiniteDifferenceStep, ...
    "GradientTolerance", options.GradientTolerance);

logResidual = log(max(path.traceResidualNorms(:), options.LogFloor));
logSolution = log(max(path.solutionNorms(:), options.LogFloor));
points = [logResidual, logSolution];

curvature = zeros(numel(alphaValues), 1);
for k = 2:(numel(alphaValues) - 1)
    left = points(k, :) - points(k - 1, :);
    right = points(k + 1, :) - points(k, :);
    chord = points(k + 1, :) - points(k - 1, :);
    denom = norm(left) * norm(right) * norm(chord);
    if denom > 0
        curvature(k) = 2 * abs(left(1) * right(2) - left(2) * right(1)) / denom;
    end
end

[maxCurvature, selectedIndex] = max(curvature);

corner = struct();
corner.kind = "fem_bem_lcurve_corner";
corner.policy = "readable_lcurve_corner_from_trace_regularization_path";
corner.path = path;
corner.alphas = alphaValues;
corner.selectedIndex = selectedIndex;
corner.selectedAlpha = alphaValues(selectedIndex);
corner.selectedTraceResidualNorm = path.traceResidualNorms(selectedIndex);
corner.selectedSolutionNorm = path.solutionNorms(selectedIndex);
corner.curvature = curvature;
corner.maxCurvature = maxCurvature;
corner.checks = struct( ...
    "pathOk", path.status == "ok", ...
    "interiorSelected", selectedIndex > 1 && selectedIndex < numel(alphaValues), ...
    "finiteCurvature", all(isfinite(curvature)), ...
    "positiveCornerCurvature", maxCurvature > 0);
corner.status = string(ifelse(all(struct2array(corner.checks)), "ok", "needs_attention"));
end


function value = ifelse(condition, a, b)
if condition
    value = a;
else
    value = b;
end
end
