function result = educationalFemBemTraceLeastSquares(model, target, options)
%EDUCATIONALFEMBEMTRACELEASTSQUARES Fit boundary data through a FEM/BEM trace.
%
% This small teaching gate connects the readable FEM/BEM mesh view to
% optimization:
%
%   phi(u) = 0.5 * (T*u - g)' * M * (T*u - g) + 0.5 * alpha * (u' * u)
%
% where T maps volume H1 P1 unknowns to boundary scalar BEM unknowns and M is
% the boundary P1 surface mass matrix.  It intentionally keeps the objective,
% gradient, normal equations, and finite-difference check visible.

arguments
    model
    target double
    options.Tikhonov (1,1) double {mustBeNonnegative} = 0
    options.Initial double = []
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.GradientTolerance (1,1) double {mustBePositive} = 1e-7
end

if ~isfield(model, "operators") || ~isfield(model.operators, "trace")
    model = assembleFirstOrderFemBemTrace(model);
end

T = model.operators.trace.matrix;
M = model.operators.bem.surfaceMass;
g = target(:);
alpha = options.Tikhonov;

nFem = size(T, 2);
if size(T, 1) ~= numel(g)
    error("educationalFemBemTraceLeastSquares:target", ...
        "Target length must match the number of scalar BEM trace nodes.");
end

if isempty(options.Initial)
    u0 = zeros(nFem, 1);
else
    u0 = options.Initial(:);
end
if numel(u0) ~= nFem
    error("educationalFemBemTraceLeastSquares:initial", ...
        "Initial must have one entry per volume H1 unknown.");
end

A = T.' * M * T + alpha * speye(nFem);
rhs = T.' * M * g;
Afull = full(A);
rankNormal = rank(Afull);
if alpha == 0 && rankNormal < nFem
    u = pinv(Afull) * rhs;
    solver = "minimum_norm_pinv_rank_deficient";
else
    u = A \ rhs;
    solver = "backslash";
end

r = T * u - g;
grad = traceGradient(T, M, g, alpha, u);
grad0 = traceGradient(T, M, g, alpha, u0);
gradFd = finiteDifferenceGradient(T, M, g, alpha, u0, options.FiniteDifferenceStep);
maxAbsGradientError = max(abs(grad0 - gradFd));

result = struct();
result.kind = "educational_fem_bem_trace_least_squares";
result.policy = "readable_trace_optimization_gate_first_order_tri_tet";
result.objective = objectiveValue(T, M, g, alpha, u);
result.objectiveAtInitial = objectiveValue(T, M, g, alpha, u0);
result.traceResidualNorm = norm(r);
result.weightedTraceResidual = sqrt(max(r.' * M * r, 0));
result.normalEquationResidual = norm(grad);
result.u = u;
result.trace = T * u;
result.target = g;
result.initial = u0;
result.tikhonov = alpha;
result.gradientCheck = struct( ...
    "point", u0, ...
    "analytic", grad0, ...
    "finiteDifference", gradFd, ...
    "step", options.FiniteDifferenceStep, ...
    "maxAbsError", maxAbsGradientError, ...
    "tolerance", options.GradientTolerance, ...
    "passed", maxAbsGradientError <= options.GradientTolerance);
result.matrix = struct( ...
    "traceRows", size(T, 1), ...
    "femUnknowns", nFem, ...
    "traceNnz", nnz(T), ...
    "massRows", size(M, 1), ...
    "rankNormal", rankNormal, ...
    "solver", solver);
end


function value = objectiveValue(T, M, g, alpha, u)
r = T * u - g;
value = 0.5 * (r.' * M * r) + 0.5 * alpha * (u.' * u);
end


function grad = traceGradient(T, M, g, alpha, u)
grad = T.' * M * (T * u - g) + alpha * u;
end


function grad = finiteDifferenceGradient(T, M, g, alpha, u, h)
grad = zeros(size(u));
for k = 1:numel(u)
    up = u;
    um = u;
    up(k) = up(k) + h;
    um(k) = um(k) - h;
    grad(k) = (objectiveValue(T, M, g, alpha, up) - objectiveValue(T, M, g, alpha, um)) / (2 * h);
end
end
