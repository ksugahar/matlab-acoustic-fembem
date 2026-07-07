function result = drumRollConvolutionQuadrature(volFile, options)
%DRUMROLLCONVOLUTIONQUADRATURE Two-spot alternating drumhead strikes (a drum roll) by CQ.
%
%   result = drumRollConvolutionQuadrature() radiates a drum roll -- a CYLINDRICAL
%   drum whose DRUMHEAD (the top face) is struck ALTERNATELY at TWO spots --
%   with the Lubich CQ time-domain BEM.  Each beat is a localized Ricker "tap"
%   on the drumhead: spot A (at +x on the head) fires on beats 1,3,5,..., spot B
%   (at -x) on beats 2,4,6,...  The only thing that changes vs a single strike is
%   BoundaryTimeData (N x nBoundary); the CQ solver carries the alternating
%   drumhead taps to the listeners.  This showcases the solver's arbitrary
%   space-time boundary excitation on a real drum shape (a cylinder, not a sphere).
%
%   Two default listeners sit above the +x and -x sides so the roll is directional:
%   the +x listener hears the A-taps (odd beats) first and louder, the -x
%   listener hears the B-taps (even beats).  result.directional records that
%   contrast; result.beatSpot is the A/B/A/... beat labelling.
%
%   result inherits every field of volTdBemConvolutionQuadrature (pressure,
%   summary, checks, status) and adds the drum-roll excitation + directionality.

arguments
    volFile (1,1) string = "S:/MATLAB/Gypsilab/fixtures/mesh_topology/drum_cylinder.vol"
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.NumBeats (1,1) double {mustBeInteger, mustBePositive} = 3
    options.BeatInterval (1,1) double {mustBePositive} = 1.0
    options.TimeStep (1,1) double {mustBePositive} = 0.3
    options.TailTime (1,1) double {mustBePositive} = 3.0
    options.StrikeWidthTime (1,1) double {mustBePositive} = 0.35
    options.StrikeWidthSpace (1,1) double {mustBePositive} = 0.33
    options.StrikeOffset (1,1) double {mustBePositive} = 0.6
    options.ListenerRadius (1,1) double {mustBePositive} = 2.5
    options.ListenerHeight (1,1) double {mustBeNonnegative} = 0.6
    options.CqRadius double = []
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder,[1 3 7])} = 1
    options.Method (1,1) string {mustBeMember(options.Method,["BDF1","BDF2"])} = "BDF2"
    options.ObservationPoints double = []
end

mesh = VolMesh(volFile);
surface = mesh.boundary();
X = surface.vtx;                                  % nBoundary x 3
nB = size(X, 1);
c0 = options.SoundSpeed;

% --- two strike spots on the DRUMHEAD (top face) at +x and -x ----------------- %
zTop = max(X(:, 3));
headTol = 0.05 * (zTop - min(X(:, 3)) + eps);     % nodes on the top drumhead
topMask = abs(X(:, 3) - zTop) < max(headTol, 1e-6);
spotA = [ options.StrikeOffset 0 zTop];           % +x drumhead spot
spotB = [-options.StrikeOffset 0 zTop];           % -x drumhead spot
wA = exp(-sum((X - spotA).^2, 2) / (2*options.StrikeWidthSpace^2)) .* topMask;
wB = exp(-sum((X - spotB).^2, 2) / (2*options.StrikeWidthSpace^2)) .* topMask;

% --- time axis long enough to hold every beat plus a radiated tail ------------ %
dt = options.TimeStep;
finalTime = options.NumBeats*options.BeatInterval + options.TailTime;
N = max(8, ceil(finalTime / dt));
t = (0:N-1).' * dt;

% --- alternating localized Ricker taps: A on odd beats, B on even beats ------- %
boundaryData = zeros(N, nB);
beatSpot = strings(options.NumBeats, 1);
beatTime = zeros(options.NumBeats, 1);
for b = 1:options.NumBeats
    tb = b * options.BeatInterval;                % beat time
    tap = rickerTap(t, tb, options.StrikeWidthTime);
    if mod(b, 2) == 1
        boundaryData = boundaryData + tap .* wA.';
        beatSpot(b) = "A";
    else
        boundaryData = boundaryData + tap .* wB.';
        beatSpot(b) = "B";
    end
    beatTime(b) = tb;
end

% --- default listeners above the +x (A) and -x (B) sides of the drum ---------- %
obs = options.ObservationPoints;
if isempty(obs)
    obs = [ options.ListenerRadius 0 zTop+options.ListenerHeight;
           -options.ListenerRadius 0 zTop+options.ListenerHeight];
end

% --- radiate the drum roll with the (mass-consistent) CQ solver --------------- %
cq = volTdBemConvolutionQuadrature(volFile, ...
    NumTime=N, TimeStep=dt, SoundSpeed=c0, Method=options.Method, ...
    QuadratureOrder=options.QuadratureOrder, CqRadius=options.CqRadius, ...
    BoundaryTimeData=boundaryData, ObservationPoints=obs);

result = cq;
result.kind = "drum_roll_two_spot_alternating_cq_time_response";
result.strikeSpotA = spotA;   result.strikeSpotB = spotB;
result.strikeWeightsA = wA;   result.strikeWeightsB = wB;
result.beatTime = beatTime;   result.beatSpot = beatSpot;
result.listeners = obs;
% mean excitation at each spot's node cluster, for a clean two-line plot
result.excitationA = boundaryData * wA / sum(wA);
result.excitationB = boundaryData * wB / sum(wB);

% --- directionality: split the odd- (A) and even- (B) beat energy ------------- %
oddBeats  = beatTime(1:2:end);
evenBeats = beatTime(2:2:end);
eA_listenerA = beatWindowEnergy(t, obs, result.pressure, 1, oddBeats,  c0, spotA, options.BeatInterval/2);
eB_listenerA = beatWindowEnergy(t, obs, result.pressure, 1, evenBeats, c0, spotB, options.BeatInterval/2);
eA_listenerB = beatWindowEnergy(t, obs, result.pressure, 2, oddBeats,  c0, spotA, options.BeatInterval/2);
eB_listenerB = beatWindowEnergy(t, obs, result.pressure, 2, evenBeats, c0, spotB, options.BeatInterval/2);
% The robust directional claim is a CROSS-LISTENER one: the A-taps are loudest
% on the A (+x) side.  (A same-listener A-vs-B split is fragile here because the
% odd/even beat counts differ -- 2 A-taps vs 1 B-tap -- so total A-energy leaks
% into every B-window; the cross-listener A-side comparison is clean.)
result.directional = struct( ...
    "listenerA_side", "+x (spot A)", ...
    "listenerB_side", "-x (spot B)", ...
    "A_energy_at_listenerA", eA_listenerA, ...
    "B_energy_at_listenerA", eB_listenerA, ...
    "A_energy_at_listenerB", eA_listenerB, ...
    "B_energy_at_listenerB", eB_listenerB, ...
    "A_louder_on_A_side", eA_listenerA > eA_listenerB, ...
    "A_side_energy_ratio", eA_listenerA / max(eA_listenerB, eps));
result.checks.alternating_two_spot = numel(unique(beatSpot)) == 2 && beatSpot(1) == "A";
result.checks.directional_drum_roll = result.directional.A_louder_on_A_side;
result.checks.drumhead_strikes_on_top_face = ...
    abs(spotA(3) - zTop) < 1e-9 && abs(spotB(3) - zTop) < 1e-9 && nnz(topMask) > 0;
end


function e = beatWindowEnergy(t, obs, pressure, listener, beatTimes, c0, spot, halfWindow)
%BEATWINDOWENERGY Energy at a listener within +/-halfWindow of each beat's arrival.
travel = norm(obs(listener,:) - spot) / c0;
sig = real(pressure(:, listener));
e = 0;
for j = 1:numel(beatTimes)
    tArr = beatTimes(j) + travel;
    win = abs(t - tArr) <= halfWindow;
    e = e + sum(sig(win).^2);
end
end


function tap = rickerTap(t, centerTime, width)
x = (t - centerTime) / width;
tap = (1 - 2*x.^2) .* exp(-x.^2);
tap(abs(tap) < 1e-14) = 0;
end
