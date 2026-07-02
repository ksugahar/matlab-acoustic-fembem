function path = femBemTikhonovPath(model, target, alphas, options)
%femBemTikhonovPath Trace-fit regularization path for FEM/BEM education.
%
% This keeps the Tikhonov trade-off visible:
%
%   phi_alpha(u) = 0.5*(T*u-g)'*M*(T*u-g) + 0.5*alpha*(u'*u).
%
% Larger alpha should damp the FEM unknown norm while usually increasing the
% boundary trace residual.  The function is intentionally a thin readable loop
% over femBemTraceLeastSquares, not a black-box optimizer.

arguments
    model
    target double
    alphas double {mustBeNonnegative}
    options.Initial double = []
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.GradientTolerance (1,1) double {mustBePositive} = 1e-7
end

alphaValues = alphas(:);
if isempty(alphaValues)
    error("femBemTikhonovPath:alphas", ...
        "At least one Tikhonov weight is required.");
end
if any(diff(alphaValues) < 0)
    error("femBemTikhonovPath:alphas", ...
        "Tikhonov weights must be sorted in nondecreasing order.");
end

rows = repmat(emptyRow(), numel(alphaValues), 1);
for k = 1:numel(alphaValues)
    fit = femBemTraceLeastSquares(model, target, ...
        "Tikhonov", alphaValues(k), ...
        "Initial", options.Initial, ...
        "FiniteDifferenceStep", options.FiniteDifferenceStep, ...
        "GradientTolerance", options.GradientTolerance);

    rows(k).alpha = alphaValues(k);
    rows(k).objective = fit.objective;
    rows(k).traceResidualNorm = fit.traceResidualNorm;
    rows(k).weightedTraceResidual = fit.weightedTraceResidual;
    rows(k).solutionNorm = norm(fit.u);
    rows(k).normalEquationResidual = fit.normalEquationResidual;
    rows(k).gradientCheckMaxAbsError = fit.gradientCheck.maxAbsError;
    rows(k).gradientCheckPassed = fit.gradientCheck.passed;
    rows(k).solver = fit.matrix.solver;
end

solutionNorms = [rows.solutionNorm].';
traceResiduals = [rows.traceResidualNorm].';
weightedResiduals = [rows.weightedTraceResidual].';

tol = 1e-10;
path = struct();
path.kind = "fem_bem_tikhonov_path";
path.policy = "readable_trace_regularization_path_first_order_tri_tet";
path.alphas = alphaValues;
path.rows = rows;
path.solutionNorms = solutionNorms;
path.traceResidualNorms = traceResiduals;
path.weightedTraceResiduals = weightedResiduals;
path.checks = struct( ...
    "alphasSorted", all(diff(alphaValues) >= 0), ...
    "solutionNormNonincreasing", all(diff(solutionNorms) <= tol), ...
    "traceResidualNondecreasing", all(diff(traceResiduals) >= -tol), ...
    "weightedResidualNondecreasing", all(diff(weightedResiduals) >= -tol), ...
    "allGradientChecksPassed", all([rows.gradientCheckPassed]));
path.status = string(ifelse(all(struct2array(path.checks)), "ok", "needs_attention"));
end


function row = emptyRow()
row = struct( ...
    "alpha", NaN, ...
    "objective", NaN, ...
    "traceResidualNorm", NaN, ...
    "weightedTraceResidual", NaN, ...
    "solutionNorm", NaN, ...
    "normalEquationResidual", NaN, ...
    "gradientCheckMaxAbsError", NaN, ...
    "gradientCheckPassed", false, ...
    "solver", "");
end


function value = ifelse(condition, a, b)
if condition
    value = a;
else
    value = b;
end
end
