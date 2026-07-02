function result = touchstoneDesignFrequencyGrid(frequencies, designFrequency, options)
%touchstoneDesignFrequencyGrid Check sweep rows around a design frequency.
%
% A Touchstone row can pass row-level passivity, reciprocity, Z0, and S11
% match gates while the surrounding frequency sweep is still too sparse for
% interpolation.  This readable MATLAB helper records the nearest row plus
% the lower/upper bracket rows before a sweep is reused in an optimization,
% equivalent-circuit, group-delay, or solver-ready notebook.

arguments
    frequencies (:,:) double {mustBeNonnegative}
    designFrequency (1,1) double {mustBePositive}
    options.FrequencyUnit (1,1) string = "GHz"
    options.DesignFrequencyUnit (1,1) string = ""
    options.MaxRelativeSpacing (1,1) double {mustBeNonnegative} = 0.05
    options.RequireBracket (1,1) logical = true
end

if isempty(frequencies)
    error("touchstoneDesignFrequencyGrid:empty", ...
        "frequencies must contain at least one row");
end
if ~isvector(frequencies)
    error("touchstoneDesignFrequencyGrid:vector", ...
        "frequencies must be a vector");
end

frequencyScale = frequencyUnitScale(options.FrequencyUnit);
if strlength(options.DesignFrequencyUnit) == 0
    designScale = frequencyScale;
    designUnit = options.FrequencyUnit;
else
    designScale = frequencyUnitScale(options.DesignFrequencyUnit);
    designUnit = options.DesignFrequencyUnit;
end

frequencyHz = frequencies(:) * frequencyScale;
designFrequencyHz = designFrequency * designScale;
if any(diff(frequencyHz) <= 0)
    error("touchstoneDesignFrequencyGrid:monotonic", ...
        "frequency grid must be strictly increasing");
end

[nearestAbsErrorHz, nearestIndex] = min(abs(frequencyHz - designFrequencyHz));
exactTolerance = max(1e-12, abs(designFrequencyHz) * 1e-12);
exactIndex = find(abs(frequencyHz - designFrequencyHz) <= exactTolerance, 1, "first");
if ~isempty(exactIndex)
    lowerIndex = exactIndex;
    upperIndex = exactIndex;
else
    lowerIndex = find(frequencyHz < designFrequencyHz, 1, "last");
    upperIndex = find(frequencyHz > designFrequencyHz, 1, "first");
end

bracketed = ~isempty(lowerIndex) && ~isempty(upperIndex);
if bracketed
    bracketGapHz = frequencyHz(upperIndex) - frequencyHz(lowerIndex);
    bracketGapRel = bracketGapHz / designFrequencyHz;
else
    bracketGapHz = NaN;
    bracketGapRel = NaN;
    lowerIndex = NaN;
    upperIndex = NaN;
end

spacingOk = bracketed && bracketGapRel <= options.MaxRelativeSpacing;
if ~options.RequireBracket && ~bracketed
    spacingOk = true;
end

checks = struct( ...
    "frequencyGridStrictlyIncreasing", true, ...
    "designFrequencyBracketed", bracketed || ~options.RequireBracket, ...
    "designSpacingOk", spacingOk, ...
    "nearestRowRecorded", ~isnan(nearestIndex));

issues = strings(0, 1);
if ~checks.designFrequencyBracketed
    issues(end + 1, 1) = "design frequency is outside the exported sweep";
end
if ~checks.designSpacingOk
    issues(end + 1, 1) = "frequency rows around the design point are too sparse";
end

result = struct();
result.kind = "touchstone_design_frequency_grid";
result.policy = "readable_touchstone_design_frequency_bracket_gate";
result.frequencyUnit = options.FrequencyUnit;
result.designFrequencyUnit = designUnit;
result.frequencyHz = frequencyHz;
result.designFrequencyHz = designFrequencyHz;
result.nRows = numel(frequencyHz);
result.nearestIndex = nearestIndex;
result.nearestFrequencyHz = frequencyHz(nearestIndex);
result.nearestAbsErrorHz = nearestAbsErrorHz;
result.nearestRelError = nearestAbsErrorHz / designFrequencyHz;
result.lowerIndex = lowerIndex;
result.upperIndex = upperIndex;
if bracketed
    result.lowerFrequencyHz = frequencyHz(lowerIndex);
    result.upperFrequencyHz = frequencyHz(upperIndex);
else
    result.lowerFrequencyHz = NaN;
    result.upperFrequencyHz = NaN;
end
result.bracketGapHz = bracketGapHz;
result.bracketGapRel = bracketGapRel;
result.maxRelativeSpacing = options.MaxRelativeSpacing;
result.requireBracket = options.RequireBracket;
result.indexBase = 1;
result.checks = checks;
result.issues = issues;
result.notes = [
    "MATLAB indices are 1-based; compare carefully with Python/radia artifacts"
    "the nearest design row is not enough unless the sweep also brackets the design frequency"
    "run this before solver-ready row preflight when interpolation or fitting will use the sweep"
];

result.status = "ok";
checkValues = struct2cell(checks);
if ~all([checkValues{:}])
    result.status = "needs_attention";
end
end


function scale = frequencyUnitScale(unit)
switch lower(strtrim(unit))
    case "hz"
        scale = 1;
    case "khz"
        scale = 1e3;
    case "mhz"
        scale = 1e6;
    case "ghz"
        scale = 1e9;
    otherwise
        error("touchstoneDesignFrequencyGrid:frequencyUnit", ...
            "Unsupported frequency unit: %s", unit);
end
end
