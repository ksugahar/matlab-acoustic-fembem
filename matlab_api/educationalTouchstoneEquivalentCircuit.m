function result = educationalTouchstoneEquivalentCircuit(s11, s21, options)
%EDUCATIONALTOUCHSTONEEQUIVALENTCIRCUIT Readable S -> Y/Z circuit extraction.
%
% This is a small MATLAB teaching companion for CST/Touchstone and
% ngsolve.bem port-mode rows.  It keeps the reference impedance visible:
%
%   Y = (I - S) * inv(I + S) / Z0
%   Z = Z0 * (I + S) * inv(I - S)
%
% Then it exposes the equivalent pi admittances and T impedances.  The same
% normalized S row produces different circuit values when Z0 changes.

arguments
    s11 (1,1) double
    s21 (1,1) double
    options.S12 (1,1) double = NaN
    options.S22 (1,1) double = NaN
    options.Z0 (1,1) double {mustBePositive} = 50
    options.ComparisonZ0 (1,1) double = NaN
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

s12 = options.S12;
s22 = options.S22;
if isnan(s12)
    s12 = s21;
end
if isnan(s22)
    s22 = s11;
end

primary = convertRow(s11, s21, s12, s22, options.Z0, options.Tolerance);

result = primary;
result.kind = "educational_touchstone_equivalent_circuit";
result.policy = "readable_touchstone_z0_equivalent_circuit_gate";
result.notes = [
    "data_format must be normalized before this MATLAB helper is called"
    "reference impedance Z0 is part of the physical circuit extraction contract"
    "for a matched isolated row S=0, y_shunt=1/Z0 and z11=Z0"
];

if ~isnan(options.ComparisonZ0)
    if options.ComparisonZ0 <= 0
        error("educationalTouchstoneEquivalentCircuit:comparisonZ0", ...
            "ComparisonZ0 must be positive when supplied.");
    end
    comparison = convertRow(s11, s21, s12, s22, options.ComparisonZ0, options.Tolerance);
    result.comparison = comparison;
    result.checks.equivalentValuesDependOnZ0 = ...
        abs(result.pi.yShunt1 - comparison.pi.yShunt1) > options.Tolerance || ...
        abs(result.z.z11 - comparison.z.z11) > options.Tolerance;
else
    result.comparison = struct();
    result.checks.equivalentValuesDependOnZ0 = true;
end

result.status = "ok";
checkNames = fieldnames(result.checks);
for k = 1:numel(checkNames)
    if ~result.checks.(checkNames{k})
        result.status = "needs_attention";
        break
    end
end
end


function row = convertRow(s11, s21, s12, s22, z0, tol)
S = [s11 s12; s21 s22];
I = eye(2);

dY = det(I + S);
if abs(dY) <= tol
    error("educationalTouchstoneEquivalentCircuit:singularY", ...
        "I+S is singular; short-circuit admittance matrix is undefined.");
end

dZ = det(I - S);
if abs(dZ) <= tol
    error("educationalTouchstoneEquivalentCircuit:singularZ", ...
        "I-S is singular; open-circuit impedance matrix is undefined.");
end

Y = ((I - S) / (I + S)) / z0;
Z = z0 * ((I + S) / (I - S));
singularValuesSquared = eig(S' * S);
maxSingularValueSquared = max(real(singularValuesSquared));
reciprocityError = abs(s21 - s12);

row = struct();
row.z0 = z0;
row.s = struct("s11", s11, "s21", s21, "s12", s12, "s22", s22);
row.y = struct( ...
    "matrix", Y, ...
    "y11", Y(1, 1), ...
    "y12", Y(1, 2), ...
    "y21", Y(2, 1), ...
    "y22", Y(2, 2));
row.z = struct( ...
    "matrix", Z, ...
    "z11", Z(1, 1), ...
    "z12", Z(1, 2), ...
    "z21", Z(2, 1), ...
    "z22", Z(2, 2));
row.pi = struct( ...
    "yShunt1", Y(1, 1) + Y(1, 2), ...
    "yShunt2", Y(2, 2) + Y(1, 2), ...
    "ySeries", -Y(1, 2));
row.t = struct( ...
    "zSeries1", Z(1, 1) - Z(1, 2), ...
    "zSeries2", Z(2, 2) - Z(1, 2), ...
    "zShunt", Z(1, 2));
row.health = struct( ...
    "reciprocityError", reciprocityError, ...
    "maxSingularValueSquared", maxSingularValueSquared, ...
    "passiveMargin", 1 - maxSingularValueSquared);
row.checks = struct( ...
    "referenceImpedancePositive", z0 > 0, ...
    "sparameterReciprocityOk", reciprocityError <= tol, ...
    "sparameterPassivityOk", maxSingularValueSquared <= 1 + tol, ...
    "yMatrixDefined", true, ...
    "zMatrixDefined", true, ...
    "yReciprocal", abs(Y(1, 2) - Y(2, 1)) <= tol, ...
    "zReciprocal", abs(Z(1, 2) - Z(2, 1)) <= tol);
end
