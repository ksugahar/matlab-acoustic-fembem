function result = educationalBoxConstrainedLeastSquares(A, b, lower, upper, options)
%EDUCATIONALBOXCONSTRAINEDLEASTSQUARES Readable projected-gradient gate.
%
% This is the smallest constrained optimization teaching primitive:
%
%   phi(x) = 0.5 * ||A*x - b||^2
%   x_{k+1} = project_box(x_k - alpha * grad phi(x_k))
%
% It returns the iterate, objective history, active bounds, and a KKT residual
% so the bound-sign convention is visible before larger FEM/BEM design
% examples are built.

arguments
    A double
    b double
    lower double
    upper double
    options.Initial double = []
    options.StepSize (1,1) double = NaN
    options.MaxIterations (1,1) double {mustBeInteger, mustBeNonnegative} = 200
    options.StepTolerance (1,1) double {mustBePositive} = 1e-12
    options.KktTolerance (1,1) double {mustBePositive} = 1e-10
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.GradientTolerance (1,1) double {mustBePositive} = 1e-7
end

b = b(:);
lower = lower(:);
upper = upper(:);
n = size(A, 2);
if size(A, 1) ~= numel(b)
    error("educationalBoxConstrainedLeastSquares:dimension", ...
        "A row count must match the number of entries in b.");
end
if numel(lower) ~= n || numel(upper) ~= n
    error("educationalBoxConstrainedLeastSquares:bounds", ...
        "lower and upper must have one entry per design variable.");
end
if any(lower > upper)
    error("educationalBoxConstrainedLeastSquares:bounds", ...
        "Each lower bound must be less than or equal to the upper bound.");
end

if isempty(options.Initial)
    rawInitial = zeros(n, 1);
else
    rawInitial = options.Initial(:);
end
if numel(rawInitial) ~= n
    error("educationalBoxConstrainedLeastSquares:initial", ...
        "Initial must have one entry per design variable.");
end

if isnan(options.StepSize)
    lipschitz = norm(A)^2;
    if lipschitz == 0
        stepSize = 1.0;
    else
        stepSize = 1.0 / lipschitz;
    end
else
    stepSize = options.StepSize;
end
if stepSize <= 0
    error("educationalBoxConstrainedLeastSquares:step", ...
        "StepSize must be positive.");
end

x = projectBox(rawInitial, lower, upper);
objectiveHistory = zeros(options.MaxIterations + 1, 1);
objectiveHistory(1) = objectiveValue(A, b, x);
lastStepNorm = 0.0;
iteration = 0;
for k = 1:options.MaxIterations
    grad = quadraticGradient(A, b, x);
    xNext = projectBox(x - stepSize * grad, lower, upper);
    lastStepNorm = norm(xNext - x);
    x = xNext;
    iteration = k;
    objectiveHistory(k + 1) = objectiveValue(A, b, x);
    if lastStepNorm <= options.StepTolerance
        break
    end
end
objectiveHistory = objectiveHistory(1:iteration + 1);

grad = quadraticGradient(A, b, x);
projectedGradientResidual = norm(x - projectBox(x - stepSize * grad, lower, upper)) / stepSize;
kkt = boxKktResidual(x, grad, lower, upper, options.KktTolerance);
grad0 = quadraticGradient(A, b, projectBox(rawInitial, lower, upper));
gradFd = finiteDifferenceGradient(A, b, projectBox(rawInitial, lower, upper), options.FiniteDifferenceStep);
maxAbsGradientError = max(abs(grad0 - gradFd));
objectiveDeltas = diff(objectiveHistory);

result = struct();
result.kind = "educational_box_constrained_least_squares";
result.policy = "readable_box_projected_gradient_gate_not_optuna_owned";
result.objective = objectiveHistory(end);
result.objectiveHistory = objectiveHistory;
result.objectiveMonotone = all(objectiveDeltas <= max(options.StepTolerance, 1e-14));
result.x = x;
result.gradient = grad;
result.initialRaw = rawInitial;
result.initialProjected = projectBox(rawInitial, lower, upper);
result.lower = lower;
result.upper = upper;
result.activeLower = x <= lower + options.KktTolerance;
result.activeUpper = x >= upper - options.KktTolerance;
result.stepSize = stepSize;
result.iterations = iteration;
result.lastStepNorm = lastStepNorm;
result.projectedGradientResidual = projectedGradientResidual;
result.kktResidual = kkt;
result.maxKktResidual = max(kkt);
result.gradientCheck = struct( ...
    "point", result.initialProjected, ...
    "analytic", grad0, ...
    "finiteDifference", gradFd, ...
    "step", options.FiniteDifferenceStep, ...
    "maxAbsError", maxAbsGradientError, ...
    "tolerance", options.GradientTolerance, ...
    "passed", maxAbsGradientError <= options.GradientTolerance);
result.matrix = struct( ...
    "rows", size(A, 1), ...
    "cols", n, ...
    "rank", rank(A), ...
    "conditionEstimate", cond(A));
end


function value = objectiveValue(A, b, x)
r = A * x - b;
value = 0.5 * (r.' * r);
end


function g = quadraticGradient(A, b, x)
g = A.' * (A * x - b);
end


function x = projectBox(x, lower, upper)
x = min(max(x, lower), upper);
end


function r = boxKktResidual(x, grad, lower, upper, tol)
r = zeros(size(x));
for k = 1:numel(x)
    if x(k) <= lower(k) + tol
        r(k) = max(0.0, -grad(k));
    elseif x(k) >= upper(k) - tol
        r(k) = max(0.0, grad(k));
    else
        r(k) = abs(grad(k));
    end
end
end


function g = finiteDifferenceGradient(A, b, x, h)
g = zeros(size(x));
for k = 1:numel(x)
    xp = x;
    xm = x;
    xp(k) = xp(k) + h;
    xm(k) = xm(k) - h;
    g(k) = (objectiveValue(A, b, xp) - objectiveValue(A, b, xm)) / (2 * h);
end
end
