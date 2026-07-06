function fig = visualizeConvolutionQuadrature(resultOrVol, options)
%VISUALIZECONVOLUTIONQUADRATURE See the Lubich CQ time-domain BEM solver work.
%
%   fig = visualizeConvolutionQuadrature() runs volTdBemConvolutionQuadrature on
%   a unit-sphere fixture and draws a six-panel "X-ray" of the method, so the
%   whole convolution-quadrature idea is visible in one figure instead of hidden
%   inside a for-loop over complex frequencies.
%
%   fig = visualizeConvolutionQuadrature("model.vol") runs the solver on that
%   mesh first.  fig = visualizeConvolutionQuadrature(result) reuses a result
%   struct you already computed (the struct returned by
%   volTdBemConvolutionQuadrature -- it already carries every internal quantity
%   drawn below: cqZeta, cqLaplaceParameter, conditionNumbers, relativeResiduals,
%   boundaryData, boundaryDensity, pressure, time).
%
%   The six panels follow the solver top to bottom, each mapped to one line of
%   the algorithm:
%
%     (1) CQ contour zeta = rho * exp(-2i pi n/N)         <- the sampling circle
%     (2) Laplace nodes  s = delta(zeta)/dt               <- WHY CQ is stable:
%                                                            every node has Re(s)>0,
%                                                            so exp(-s r/c) DECAYS
%                                                            (a screened-Laplace
%                                                            kernel), never an
%                                                            oscillatory Helmholtz
%                                                            one.  A-stability is
%                                                            just "the contour
%                                                            stays in the right
%                                                            half-plane".
%     (3) cond(V(s)) and residual per node                <- N independent, well
%                                                            posed frequency solves
%     (4) source pulse g(t) and its spectrum              <- the FFT input
%     (5) surface density q(t) as a space-time map        <- the solved unknown,
%                                                            recovered by the iFFT
%                                                            with the rho^-n unscaling
%     (6) pressure impulse response p(x_obs, t)           <- the physics you can see:
%                                                            a surface pulse radiates
%                                                            a delayed impulse outward
%
%   Name-value options:
%     Result / Vol handled positionally above.
%     SavePath   (string)  export the figure to a PNG (uses exportgraphics).
%     NumTime, TimeStep, SoundSpeed, Method, QuadratureOrder  forwarded to the
%                          solver when a .vol / default run is requested.
%     Verbose    (logical) print the step-by-step narration to the console
%                          (default true) so running the file is itself a lesson.

arguments
    resultOrVol = ""
    options.SavePath (1,1) string = ""
    % Showcase defaults chosen to stay in the CQ well-recovered regime: N large
    % enough for a smooth space-time picture, but not so large that the rho^-n
    % unscaling amplifies round-off into a late-time blow-up (verified clean on
    % the unit-sphere fixture: max|p| ~ 0.4, no late-time growth).
    options.NumTime (1,1) double {mustBeInteger, mustBeGreaterThan(options.NumTime, 3)} = 24
    options.TimeStep (1,1) double {mustBePositive} = 0.2
    options.SoundSpeed (1,1) double {mustBePositive} = 1.0
    options.Method (1,1) string {mustBeMember(options.Method, ["BDF1", "BDF2"])} = "BDF2"
    options.QuadratureOrder (1,1) double {mustBeMember(options.QuadratureOrder, [1 3 7])} = 3
    options.Verbose (1,1) logical = true
end

% ------------------------------------------------------------------ %
% 0. Obtain a result struct (run the solver unless one was handed in).
% ------------------------------------------------------------------ %
if isstruct(resultOrVol)
    result = resultOrVol;
else
    volFile = string(resultOrVol);
    if strlength(volFile) == 0
        volFile = defaultSphereFixture();
    end
    result = volTdBemConvolutionQuadrature(volFile, ...
        NumTime=options.NumTime, TimeStep=options.TimeStep, ...
        SoundSpeed=options.SoundSpeed, Method=options.Method, ...
        QuadratureOrder=options.QuadratureOrder);
end

N      = numel(result.time);
zeta   = result.cqZeta(:);
s      = result.cqLaplaceParameter(:);
rho    = result.cqRadius;
t      = result.time(:);
gNode  = result.boundaryData(:, 1);                 % source pulse at node 1
density = result.boundaryDensity;                   % N x nBoundary, real
pressure = result.pressure;                         % N x nObs, real
nodeIdx = (0:N-1).';

say(options.Verbose, "");
say(options.Verbose, "Convolution-quadrature X-ray (%d time steps, %s, c = %g)", ...
    N, result.method, options.SoundSpeed);
say(options.Verbose, "  status = %s", result.status);

% ------------------------------------------------------------------ %
% Figure scaffold (batch-safe: invisible when we are only saving a PNG).
% ------------------------------------------------------------------ %
visible = "on";
if strlength(options.SavePath) > 0 && ~usejava("desktop")
    visible = "off";
end
fig = figure("Color", "w", "Position", [80 80 1280 760], "Visible", visible);
tl = tiledlayout(fig, 2, 3, "TileSpacing", "compact", "Padding", "compact");
title(tl, "Lubich convolution quadrature, step by step", ...
    "FontWeight", "bold", "FontSize", 13);

% ------------------------------------------------------------------ %
% (1) The sampling contour zeta on the rho-circle.
% ------------------------------------------------------------------ %
say(options.Verbose, "(1) sample zeta on a circle of radius rho = %.4g < 1", rho);
ax1 = nexttile(tl, 1);
drawUnitCircle(ax1, 1.0, "k:");
drawUnitCircle(ax1, rho, [0.6 0.6 0.6]);
scatter(ax1, real(zeta), imag(zeta), 36, nodeIdx, "filled");
axis(ax1, "equal"); grid(ax1, "on"); styleAxes(ax1);
xlabel(ax1, "Re \zeta"); ylabel(ax1, "Im \zeta");
title(ax1, {"(1) CQ contour", sprintf("\\zeta = \\rho e^{-2\\pi i n/N},  \\rho = %.3g", rho)});

% ------------------------------------------------------------------ %
% (2) The Laplace-domain nodes s = delta(zeta)/dt -- the whole point.
% ------------------------------------------------------------------ %
say(options.Verbose, "(2) map to s = delta(zeta)/dt   ->   min Re(s) = %.4g (must be > 0)", min(real(s)));
ax2 = nexttile(tl, 2);
xlims = [min(real(s)) max(real(s))]; ylims = [min(imag(s)) max(imag(s))];
padx = 0.08*range(xlims) + eps; pady = 0.08*range(ylims) + eps;
shadeRightHalfPlane(ax2, [xlims(1)-padx xlims(2)+padx], [ylims(1)-pady ylims(2)+pady]);
scatter(ax2, real(s), imag(s), 40, nodeIdx, "filled");
xline(ax2, 0, "k-", "LineWidth", 1);
grid(ax2, "on"); styleAxes(ax2);
xlabel(ax2, "Re s"); ylabel(ax2, "Im s");
title(ax2, {"(2) Laplace nodes  s = \delta(\zeta)/dt", "all Re s > 0  \rightarrow  exp(-s r/c) decays (A-stable)"});
cb2 = colorbar(ax2); ylabel(cb2, "frequency node n"); set(cb2, "Color", "k");

% ------------------------------------------------------------------ %
% (3) Conditioning + residual of the N frequency solves.
% ------------------------------------------------------------------ %
say(options.Verbose, "(3) N=%d frequency solves: max cond(V) = %.3g, max residual = %.2e", ...
    N, max(result.conditionNumbers), max(result.relativeResiduals));
ax3 = nexttile(tl, 3);
yyaxis(ax3, "left");
semilogy(ax3, nodeIdx, result.conditionNumbers, "-o", "LineWidth", 1.2, "MarkerSize", 4);
ylabel(ax3, "cond V(s)");
yyaxis(ax3, "right");
semilogy(ax3, nodeIdx, max(result.relativeResiduals, 1e-18), "-s", "LineWidth", 1.2, "MarkerSize", 4);
ylabel(ax3, "relative residual");
grid(ax3, "on"); styleAxes(ax3); xlabel(ax3, "frequency node n");
title(ax3, {"(3) conditioning + residual", "each node = one well-posed screened-Laplace solve"});

% ------------------------------------------------------------------ %
% (4) Source pulse in time + its spectrum (the FFT input).
% ------------------------------------------------------------------ %
say(options.Verbose, "(4) source pulse g(t): the RHS fed into the FFT");
ax4 = nexttile(tl, 4);
plot(ax4, t, gNode, "-o", "LineWidth", 1.4, "MarkerSize", 4, "Color", [0.10 0.35 0.70]);
grid(ax4, "on"); styleAxes(ax4); xlabel(ax4, "time t"); ylabel(ax4, "g(t) at a node");
title(ax4, "(4) source pulse  g(t)  (FFT input)");

% ------------------------------------------------------------------ %
% (5) Surface density q(t): the recovered unknown, space-time.
% ------------------------------------------------------------------ %
say(options.Verbose, "(5) surface density q(t) recovered by iFFT (rho^-n unscaling)");
ax5 = nexttile(tl, 5);
imagesc(ax5, 1:size(density, 2), t, density);
set(ax5, "YDir", "normal"); styleAxes(ax5); colormap(ax5, "turbo");
% Robust color limit: keep the physical pulse readable, do not let a few
% late-time round-off outliers (the rho^-n amplification tail) stretch the scale.
clip = prctile(abs(density(:)), 97);
if clip > 0
    clim(ax5, [-clip clip]);
end
xlabel(ax5, "boundary node"); ylabel(ax5, "time t");
title(ax5, "(5) solved surface density  q(t)");
cb5 = colorbar(ax5); ylabel(cb5, "density"); set(cb5, "Color", "k");

% ------------------------------------------------------------------ %
% (6) Pressure impulse response -- the physics you can see.
% ------------------------------------------------------------------ %
say(options.Verbose, "(6) pressure impulse response p(x_obs, t) at %d observers", size(pressure, 2));
ax6 = nexttile(tl, 6);
hold(ax6, "on");
co = lines(max(size(pressure, 2), 1));
for j = 1:size(pressure, 2)
    plot(ax6, t, pressure(:, j), "-", "LineWidth", 1.4, "Color", co(j, :), ...
        "DisplayName", sprintf("obs %d", j));
end
hold(ax6, "off");
grid(ax6, "on"); styleAxes(ax6); xlabel(ax6, "time t"); ylabel(ax6, "pressure p(x,t)");
title(ax6, "(6) radiated impulse response  p(x_{obs}, t)");
legend(ax6, "Location", "best", "TextColor", "k");

% ------------------------------------------------------------------ %
% Save if requested.
% ------------------------------------------------------------------ %
if strlength(options.SavePath) > 0
    exportgraphics(fig, options.SavePath, "Resolution", 150);
    say(options.Verbose, "saved figure -> %s", options.SavePath);
end
end


% ====================================================================== %
% local helpers (kept tiny and inline so the file reads top to bottom)
% ====================================================================== %
function say(verbose, fmt, varargin)
if verbose
    fprintf(fmt + "\n", varargin{:});
end
end


function drawUnitCircle(ax, r, spec)
th = linspace(0, 2*pi, 256);
hold(ax, "on");
if ischar(spec) || isstring(spec)
    plot(ax, r*cos(th), r*sin(th), spec, "LineWidth", 1);
else
    plot(ax, r*cos(th), r*sin(th), "-", "Color", spec, "LineWidth", 1);
end
end


function shadeRightHalfPlane(ax, xr, yr)
% Shade the stable Re(s) > 0 half-plane (light green) behind the s nodes.
x0 = max(xr(1), 0);
hold(ax, "on");
patch(ax, [x0 xr(2) xr(2) x0], [yr(1) yr(1) yr(2) yr(2)], ...
    [0.85 0.95 0.85], "EdgeColor", "none");
xlim(ax, xr); ylim(ax, yr);
end


function styleAxes(ax)
set(ax, "Color", "w", "XColor", "k", "YColor", "k", "FontSize", 11, "Box", "on");
end


function volFile = defaultSphereFixture()
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
volFile = string(fullfile(repoRoot, "fixtures", "mesh_topology", "unit_sphere_coarse.vol"));
end
