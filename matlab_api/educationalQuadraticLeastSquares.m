function result = educationalQuadraticLeastSquares(A, b, options)
%EDUCATIONALQUADRATICLEASTSQUARES Readable quadratic least-squares gate.
%
% This is a small optimization teaching primitive:
%
%   phi(x) = 0.5 * ||A*x - b||^2
%   grad phi(x) = A' * (A*x - b)
%
% It returns the closed-form least-squares solution together with an analytic
% vs central-finite-difference gradient check at the chosen initial point.

arguments
    A double
    b double
    options.Initial double = []
    options.FiniteDifferenceStep (1,1) double {mustBePositive} = 1e-6
    options.GradientTolerance (1,1) double {mustBePositive} = 1e-7
end

b = b(:);
n = size(A, 2);
if size(A, 1) ~= numel(b)
    error("educationalQuadraticLeastSquares:dimension", ...
        "A row count must match the number of entries in b.");
end

if isempty(options.Initial)
    x0 = zeros(n, 1);
else
    x0 = options.Initial(:);
end
if numel(x0) ~= n
    error("educationalQuadraticLeastSquares:initial", ...
        "Initial must have one entry per design variable.");
end

x = A \ b;
r = A * x - b;
g = quadraticGradient(A, b, x);
g0 = quadraticGradient(A, b, x0);
gfd = finiteDifferenceGradient(A, b, x0, options.FiniteDifferenceStep);
maxAbsGradientError = max(abs(g0 - gfd));

result = struct();
result.kind = "educational_quadratic_least_squares";
result.policy = "readable_matlab_optimization_gate_not_optuna_owned";
result.objective = objectiveValue(A, b, x);
result.objectiveAtInitial = objectiveValue(A, b, x0);
result.residualNorm = norm(r);
result.normalEquationResidual = norm(g);
result.x = x;
result.initial = x0;
result.gradientCheck = struct( ...
    "point", x0, ...
    "analytic", g0, ...
    "finiteDifference", gfd, ...
    "step", options.FiniteDifferenceStep, ...
    "maxAbsError", maxAbsGradientError, ...
    "tolerance", options.GradientTolerance, ...
    "passed", maxAbsGradientError <= options.GradientTolerance);
result.matrix = struct( ...
    "rows", size(A, 1), ...
    "cols", size(A, 2), ...
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
