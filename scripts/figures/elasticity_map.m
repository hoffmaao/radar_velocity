function OUT = elasticity_map(opts)
%ELASTICITY_MAP E* along the line, where the data can actually resolve it.
%
%   Slides a 1 km patch along each usable line and refits the modulus
%   INSIDE the patch only (vdef.invertElasticModulus E_patch mode), with
%   everything else held at the line's global joint fit: E_ref outside the
%   patch, the clamp fixed at the global x0, both observables in the
%   misfit. The result is a resolution-honest local E*(x): per-patch
%   profile-likelihood intervals, and NaN bounds wherever the data cannot
%   know the answer.
%
%   WHAT "APPROPRIATE" MEANS HERE, made visible rather than argued:
%     - Flexure is NONLOCAL. A patch narrower than the ~1.2 km flexural
%       length reports a smoothed average, so neighbouring patches are
%       correlated and the profile has ~3 independent elements at most.
%     - The constraint DIES SEAWARD. Where the shelf floats flat, bending
%       moment and curvature both vanish and D = M/kappa is 0/0; patches
%       there come back unconstrained (open markers, no whiskers), which
%       is the correct answer and the reason a full-line "map" would
%       overstate what the survey knows.
%     - A CONSTANT-TRUTH CONTROL runs beside the data: synthetic
%       observations from the line's own global fit with the line's own
%       sigmas, through the identical patch pipeline. Structure in the
%       real profile is only worth discussing where it exceeds what the
%       control invents from noise. The control bounds the INVERSION's
%       false structure, not the data's systematics - with chi2_strain
%       above 1 the real scatter is expected to exceed it.
%
%   Fits come from scripts/diagnostics/elastic_modulus.m via opts.fit
%   (pass the saved OUT to skip the refit), same as the other figures.
%
%   opts, all optional:
%     .fit         OUT struct from elastic_modulus
%     .out_dir     where the png goes
%     .patch_w     patch width [m] (default 1000)
%     .mp_dir, .products, .block   passed through when refitting
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/figures'); elasticity_map"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root);
addpath(fullfile(root,'scripts','diagnostics'));

if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures_flexure';
end
if ~isfield(opts,'patch_w') || isempty(opts.patch_w), opts.patch_w = 1000; end
if ~isfield(opts,'fit'), opts.fit = []; end

if isempty(opts.fit)
  dopts = struct();
  for f = {'mp_dir','products','block'}
    if isfield(opts, f{1}), dopts.(f{1}) = opts.(f{1}); end
  end
  D = elastic_modulus(dopts);
else
  D = opts.fit;
end
L = D.lines; R = D.fits; CASE = 1;
nL = numel(L);

% Same gate as elasticity_results: bounded AND the beam describes the
% line's own flexure data.
con = false(1,nL);
for i = 1:nL
  Ri = R{CASE,i};
  if isempty(Ri), continue; end
  x2a = Ri.chi2red;
  if Ri.has_strain && isfinite(Ri.chi2_shape), x2a = Ri.chi2_shape; end
  con(i) = isfinite(Ri.E_hi) && isfinite(x2a) && x2a < 3;
end

rand('seed', 7); randn('seed', 7);   %#ok<RAND> % control is reproducible

P = struct('name',{},'xc',{},'E',{},'lo',{},'hi',{},'Ec',{},'glob',{}, ...
           'lat',{},'lon',{},'x0',{});
for i = 1:nL
  if ~con(i), continue; end
  Rg = R{CASE,i};
  ok = isfinite(L(i).a) & isfinite(L(i).a_std) & L(i).a_std > 0;
  xa = L(i).x_sea(ok); ya = L(i).a(ok); sa = L(i).a_std(ok);

  base = struct('h', Rg.h_spec, 'sigma', sa, 'avg_width', 500, ...
    'E_grid', logspace(log10(Rg.E/12), log10(Rg.E*12), 31), ...
    'x0_grid', Rg.x0);
  hasS = Rg.has_strain;
  if hasS
    base.strain = struct('x', L(i).s_x, 'y', L(i).s_y, ...
      'sigma', L(i).s_sig, 'ref_depth', 100, 'avg_width', 500);
  end

  % Control observations: the global fit's own predictions plus the
  % line's own noise, so control and data differ ONLY in whether the ice
  % put structure there.
  yc = Rg.w_model + sa.*randn(size(sa));
  if hasS
    sc = Rg.strain_model + L(i).s_sig.*randn(size(L(i).s_sig));
  end

  xc = xa(:).';
  Ei = nan(size(xc)); loi = Ei; hii = Ei; Eci = Ei;
  for p = 1:numel(xc)
    po = base;
    po.E_patch = struct('lo', xc(p)-opts.patch_w/2, ...
                        'hi', xc(p)+opts.patch_w/2, 'E_ref', Rg.E);
    Rp = vdef.invertElasticModulus(xa, ya, po);
    Ei(p) = Rp.E; loi(p) = Rp.E_lo; hii(p) = Rp.E_hi;
    pc = po; pc.sigma = sa;
    if hasS, pc.strain.y = sc; end
    Rc = vdef.invertElasticModulus(xa, yc, pc);
    Eci(p) = Rc.E;
  end
  P(end+1) = struct('name', L(i).name, 'xc', xc, 'E', Ei, ...
    'lo', loi, 'hi', hii, 'Ec', Eci, 'glob', Rg.E, ...
    'lat', L(i).lat(ok), 'lon', L(i).lon(ok), 'x0', Rg.x0); %#ok<AGROW>

  fprintf('\n=== %s: local E* (patch %.0f m, clamp fixed %.2f km, E_ref %.2f GPa) ===\n', ...
    strrep(L(i).name,'EAGER_2022_',''), opts.patch_w, Rg.x0/1e3, Rg.E/1e9);
  fprintf('%9s %8s %8s %8s %10s\n','x_sea km','E GPa','lo','hi','control');
  for p = 1:numel(xc)
    fprintf('%9.2f %8.2f %8.2f %8.2f %10.2f\n', xc(p)/1e3, ...
      Ei(p)/1e9, loi(p)/1e9, hii(p)/1e9, Eci(p)/1e9);
  end
end
assert(~isempty(P), 'no line passed the constraint gate');
OUT = P;

%% Draw
cols = [0.11 0.42 0.69; 0.89 0.47 0.10; 0.20 0.60 0.25; 0.75 0.20 0.30];
linecol = @(nm) cols(mod(find(strcmp({'EAGER_2022_GL1','EAGER_2022_GL2', ...
  'EAGER_2022_GL3','EAGER_2022_GL4'}, nm), 1)-1, size(cols,1))+1, :);
grey = [0.62 0.62 0.62];

hf = figure('Visible','off','Position',[100 100 900 640]);
set(0,'CurrentFigure',hf);

% (a) the local modulus, with unresolved patches shown as exactly that
ax1 = axes('parent',hf,'Position',[0.09 0.42 0.87 0.52]);
hold(ax1,'on'); hleg = []; lleg = {};
YL = [0.2 60];
for k = 1:numel(P)
  cc = linecol(P(k).name);
  for p = 1:numel(P(k).xc)
    xk = P(k).xc(p)/1e3;
    if isfinite(P(k).lo(p)) && isfinite(P(k).hi(p))
      plot(ax1, [xk xk], [P(k).lo(p) P(k).hi(p)]/1e9, '-', ...
        'Color', cc, 'LineWidth', 1.5);
      plot(ax1, xk, P(k).E(p)/1e9, 'o', 'Color', cc, ...
        'MarkerSize', 5, 'MarkerFaceColor', cc);
    else
      % Unconstrained: the honest rendering is an open marker with no
      % whisker, not a bar to the axis that reads as an interval.
      plot(ax1, xk, P(k).E(p)/1e9, 'o', 'Color', grey, 'MarkerSize', 5);
    end
  end
  hleg(end+1) = plot(ax1, P(k).xc/1e3, P(k).Ec/1e9, ':', ...
    'Color', cc, 'LineWidth', 1.1); %#ok<AGROW>
  lleg{end+1} = sprintf('%s control', strrep(P(k).name,'EAGER_2022_','')); %#ok<AGROW>
  plot(ax1, [min(P(k).xc) max(P(k).xc)]/1e3, [1 1]*P(k).glob/1e9, '-', ...
    'Color', [cc 0.35], 'LineWidth', 0.9);
end
set(ax1,'YScale','log'); ylim(ax1, YL); grid(ax1,'on'); box(ax1,'on');
ylabel(ax1,'Local E* (GPa)');
set(ax1,'XTickLabel',[]);
title(ax1,'(a)  E* per 1 km patch, clamp and far field held at the global fit');
legend(ax1, hleg, lleg, 'Location','NorthEast');
% Excluded lines carry no patch values to draw, but absence without trace
% is the reading the map view's grey convention exists to prevent - so the
% profile panel says so in words instead of fabricating markers.
excl = {};
for i = 1:nL
  if ~con(i) && ~isempty(R{CASE,i})
    excl{end+1} = strrep(L(i).name,'EAGER_2022_',''); %#ok<AGROW>
  end
end
if ~isempty(excl)
  text(ax1, 0.02, 0.06, sprintf('%s gated out - no patch values fitted', ...
    strjoin(excl, ', ')), 'Units','normalized', 'FontSize', 8.5, ...
    'Color', grey);
  fprintf('profile panels: %s gated out - annotated, no patch values\n', ...
    strjoin(excl, ', '));
end

% (b) the resolution profile: how much the data say, patch by patch
ax2 = axes('parent',hf,'Position',[0.09 0.09 0.87 0.27]);
hold(ax2,'on');
for k = 1:numel(P)
  cc = linecol(P(k).name);
  spand = log10(P(k).hi ./ P(k).lo);
  fin = isfinite(spand);
  plot(ax2, P(k).xc(fin)/1e3, spand(fin), 'o-', 'Color', cc, ...
    'MarkerSize', 4, 'MarkerFaceColor', cc, 'LineWidth', 1.1);
end
grid(ax2,'on'); box(ax2,'on');
xlabel(ax2,'Seaward distance (km)');
ylabel(ax2,'Interval width (decades)');
title(ax2,'(b)  Where the line resolves stiffness at all');
xl = get(ax1,'XLim'); set(ax2,'XLim', xl);

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir,'EAGER_2022_elasticity_map.png');
print(hf, out_fn, '-dpng', '-r140');
fprintf('\nWrote %s\n', out_fn);

%% Map view
% EPSG:3031 Antarctic Polar Stereographic in km - the same frame, style
% and conventions as admittance_map.m (marker SHAPE is line identity,
% marker FILL is the value, north points roughly DOWN at lon ~168 E).
% Fill is the log-ratio of local E* to the pooled global, diverging about
% zero, so the map reads as softer-than / stiffer-than the line mean and
% clips rather than stretches for the unconstrained far end. Patches
% whose interval does not close are open grey - drawn, not valued. The
% MEaSUREs grounding line is the black curve.
PAL.cat_mk = {'o','s','^','d','v'};
PAL.ink = [0 0 0]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839];
MKIDX = struct('EAGER_2022_GL1',2,'EAGER_2022_GL2',3, ...
               'EAGER_2022_GL3',4,'EAGER_2022_GL4',5);
Eref = mean([P.glob]);
fprintf('map view: diverging about the pooled global %.2f GPa\n', Eref/1e9);

ps = projcrs(3031);
h2 = figure('Visible','off','Position',[100 100 640 640],'Color','w');
set(0,'CurrentFigure',h2);
axm = axes('parent',h2,'Position',[0.11 0.10 0.80 0.86]);
hold(axm,'on');

CLIMD = 0.5;                     % +/- half a decade about the pooled mean
ndiv = 256; half = ndiv/2;
dmap = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,half)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,ndiv-half))];

hTrk = gobjects(1, numel(P));
for k = 1:numel(P)
  [xk, yk] = ps_km(ps, P(k).lon, P(k).lat);
  hTrk(k) = plot(axm, xk, yk, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
  mk = PAL.cat_mk{MKIDX.(P(k).name)};
  for p = 1:numel(P(k).E)
    if isfinite(P(k).lo(p)) && isfinite(P(k).hi(p))
      v  = log10(P(k).E(p)/Eref);       % softer = red, stiffer = blue
      ci = max(1, min(ndiv, round((v+CLIMD)/(2*CLIMD)*(ndiv-1))+1));
      % Size encodes CONSTRAINT: the interval width in decades, mapped so
      % a tight patch (0.2 dec) draws at 11 pt and a barely-bounded one
      % (1.5 dec) at 5 pt. Without this the wide-interval mid-line
      % patches - red, eye-catching, and exactly the ones the control
      % says are unresolved - read as the most confident on the map.
      spand = log10(P(k).hi(p)/P(k).lo(p));
      msz = 11 - 6*min(max((spand - 0.2)/1.3, 0), 1);
      plot(axm, xk(p), yk(p), mk, 'MarkerSize', msz, ...
        'MarkerFaceColor', dmap(ci,:), 'MarkerEdgeColor', PAL.ink_soft, ...
        'LineWidth', 0.6);
    else
      plot(axm, xk(p), yk(p), mk, 'MarkerSize', 5, ...
        'MarkerEdgeColor', [0.65 0.65 0.65], 'LineWidth', 0.8);
    end
  end
end
% The fitted clamps are deliberately NOT drawn. They sit beyond the end of
% the surveyed track, so placing them needs a straight-bearing
% extrapolation, and as unlabelled symbols they read as observations at
% that spot when they are fitted parameters. The clamp-vs-grounding-line
% comparison lives in the printed table (x0 per line) and the profile
% panels instead.

% Excluded legs are DRAWN, in grey with open markers and no values - the
% same shown-not-averaged convention as the results forest. A leg absent
% without trace reads as missing data; grey presence records that it
% exists and that the gate held it out (GL2: a(x) misfits the beam at
% chi2 3.7 on the ref_z anomaly, while its strain fits at 0.65).
for i = 1:nL
  if con(i) || isempty(R{CASE,i}), continue; end
  okx = isfinite(L(i).a) & isfinite(L(i).lat) & isfinite(L(i).lon);
  if ~any(okx), continue; end
  [xk, yk] = ps_km(ps, L(i).lon(okx), L(i).lat(okx));
  hTrk(end+1) = plot(axm, xk, yk, '-', 'Color', [0.85 0.85 0.85], ...
    'LineWidth', 0.5); %#ok<AGROW>
  mk = PAL.cat_mk{MKIDX.(L(i).name)};
  plot(axm, xk, yk, mk, 'MarkerSize', 5, ...
    'MarkerEdgeColor', [0.65 0.65 0.65], 'LineWidth', 0.8);
  fprintf('map view: %s drawn as excluded (no patch values)\n', ...
    strrep(L(i).name,'EAGER_2022_',''));
end

% ApRES sites, from the same file the other maps use. Guarded like
% admittance_map's read_apres_xy: a missing or malformed sites file is a
% printed note, not an abort after the expensive patch scan.
apres_xy = '/kucresis/scratch/hoffmana_sta/vvel/gis/eastwind_2022_2023_apres_xy.txt';
if ~exist(apres_xy,'file')
  fprintf('ApRES site coordinates not found (%s) - sites not drawn.\n', apres_xy);
else
  try
    fid = fopen(apres_xy,'r'); Csx = textscan(fid,'%f %f %s'); fclose(fid);
    for q = 1:numel(Csx{3})
      plot(axm, Csx{1}(q)/1e3, Csx{2}(q)/1e3, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor','w', 'MarkerEdgeColor', PAL.ink, 'LineWidth', 1.4);
      text(axm, Csx{1}(q)/1e3 + 0.12, Csx{2}(q)/1e3, Csx{3}{q}, ...
        'FontSize', 8, 'Color', PAL.ink, 'VerticalAlignment','middle');
    end
    fprintf('ApRES sites drawn: %s\n', strjoin(Csx{3}.', ', '));
  catch ME
    fprintf('Could not read %s (%s) - sites not drawn.\n', apres_xy, ME.message);
  end
end

axis(axm,'equal');
xl = xlim(axm); ylm = ylim(axm);
xlim(axm, xl + 0.10*diff(xl)*[-1 1]);
ylim(axm, ylm + 0.10*diff(ylm)*[-1 1]);
gl_drawn = overlay_gl(axm, '/kucresis/scratch/hoffmana_sta/vvel/gis');
if gl_drawn, fprintf('map view: grounding line drawn\n'); end
grid(axm,'on');
set(axm,'GridAlpha',0.15,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off');
xlabel(axm,'Polar stereographic x (km, EPSG:3031)','Color',PAL.ink);
ylabel(axm,'Polar stereographic y (km, EPSG:3031)','Color',PAL.ink);
colormap(axm, dmap); caxis(axm, [-CLIMD CLIMD]);
cb = colorbar(axm);
% Ticks in GPa at round values, placed at their log-ratio positions. The
% candidates are filtered against the DATA-DEPENDENT Eref so the labels
% always land inside the fixed caxis - hardcoding one tick set would leave
% a nearly bare colorbar the first time a rebuild moves the pooled mean.
cand = [0.3 0.5 0.7 1 1.5 2 3 4 5 7 10 15 20 30];
tv = cand(abs(log10(cand/(Eref/1e9))) <= CLIMD*0.999);
if numel(tv) < 3
  tv = (Eref/1e9) * 10.^linspace(-0.8*CLIMD, 0.8*CLIMD, 5);
end
set(cb, 'Ticks', log10(tv/(Eref/1e9)), ...
        'TickLabels', arrayfun(@(v) sprintf('%.3g',v), tv, 'uni', 0));
set(get(cb,'ylabel'),'string','Local E* (GPa)','Color',PAL.ink);
set(cb,'XColor',PAL.ink,'YColor',PAL.ink);

rema_tif = '/kucresis/scratch/hoffmana_sta/vvel/gis/rema/17_33_10m_v2.0_browse.tif';
if rema_underlay(axm, rema_tif)
  set(hTrk, 'Color', 'w', 'LineWidth', 0.7);
  grid(axm, 'off');
end

out_fn2 = fullfile(opts.out_dir,'EAGER_2022_elasticity_map_view.png');
print(h2, out_fn2, '-dpng', '-r140');
fprintf('Wrote %s\n', out_fn2);

end

%% ========================================================================
% The three map helpers below are copied verbatim from admittance_map.m,
% following the project's existing convention of duplicating
% rema_underlay per figure script rather than sharing a path-dependent
% helper (see also tidal_deformation.m, eastwind_survey_movie.m).

%% ========================================================================
function [xk, yk] = ps_km(ps, lon, lat)
%PS_KM Project lon/lat to EPSG:3031 Antarctic Polar Stereographic, in km.
[x, y] = projfwd(ps, lat, lon);
xk = x/1e3; yk = y/1e3;
end

%% ========================================================================
function ok = rema_underlay(ax, tif)
%REMA_UNDERLAY Draw the REMA v2 hillshade under an EPSG:3031 axes in km.
%   The axes are already in the raster's own projection, so the block
%   covering the current view is read and dropped straight in - no
%   resampling, no rotation, no reprojection. It is lifted into a light
%   gray range so the data drawn on top stays dominant, and pushed to the
%   bottom of the draw order. The block is padded 8% beyond the limits so
%   a small later limit adjustment does not expose white strips. Returns
%   false (with a message) when the tile is missing or unreadable, and
%   the figure then renders exactly as it did without imagery.
%
%   The browse tif is a COG whose overview IFDs make geotiffinfo error
%   ("multiple images ... sizes are different"), so the georeference comes
%   from georasterinfo and the pixels from imread on IFD 1.
ok = false;
if ~exist(tif,'file')
  fprintf('REMA hillshade not found (%s) - no imagery underlay.\n', tif);
  return;
end
try
  R3 = georasterinfo(tif).RasterReference;
  xl = xlim(ax); yl = ylim(ax);                     % km, EPSG:3031
  xg = xl + 0.08*diff(xl)*[-1 1];
  yg = yl + 0.08*diff(yl)*[-1 1];
  px = R3.CellExtentInWorldX;
  py = R3.CellExtentInWorldY;                       % may differ from px
  c0 = max(1, floor((xg(1)*1e3 - R3.XWorldLimits(1))/px));
  c1 = min(R3.RasterSize(2), ceil((xg(2)*1e3 - R3.XWorldLimits(1))/px));
  r0 = max(1, floor((R3.YWorldLimits(2) - yg(2)*1e3)/py));
  r1 = min(R3.RasterSize(1), ceil((R3.YWorldLimits(2) - yg(1)*1e3)/py));
  if c1 <= c0 || r1 <= r0
    fprintf('REMA tile does not cover this view - no imagery underlay.\n');
    return;
  end
  A = double(imread(tif, 'Index', 1, 'PixelRegion', {[r0 r1],[c0 c1]}));
  A(A == 0) = NaN;                                  % nodata
  vv = sort(A(isfinite(A)));
  if isempty(vv)
    fprintf('REMA tile is empty over this view - no imagery underlay.\n');
    return;
  end
  lo = vv(max(1,round(0.02*numel(vv)))); hi = vv(round(0.98*numel(vv)));
  g = min(max((A - lo)/max(hi-lo, 1), 0), 1);
  g = 0.58 + 0.40*g;                                % recessive light grays
  g(~isfinite(A)) = 1;                              % nodata as paper white
  % cell centres of the block actually read, in km; y descends with row
  ximg = (R3.XWorldLimits(1) + ([c0 c1] - 0.5)*px)/1e3;
  yimg = (R3.YWorldLimits(2) - ([r0 r1] - 0.5)*py)/1e3;
  hImg = image(ax, 'XData', ximg, 'YData', yimg, 'CData', repmat(g,[1 1 3]));
  uistack(hImg, 'bottom');
  xlim(ax, xl); ylim(ax, yl);
  ok = true;
  fprintf('REMA hillshade underlay: %d x %d px, native EPSG:3031\n', ...
    r1-r0+1, c1-c0+1);
catch ME
  fprintf('REMA underlay failed (%s) - continuing without imagery.\n', ME.message);
end
end

%% ========================================================================
function drawn = overlay_gl(ax, gis_dir)
drawn = false;
fn = fullfile(gis_dir,'GroundingLine_Antarctica_v02.shp');
if ~exist(fn,'file'), return; end
try
  S = shaperead(fn);
  X = []; Y = [];
  for k = 1:numel(S)
    X = [X; S(k).X(:); NaN]; Y = [Y; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  gx = X/1e3; gy = Y/1e3;      % the shapefile is already EPSG:3031 metres
  xl = xlim(ax); yl = ylim(ax); pad = 3;
  bad = gx < xl(1)-pad | gx > xl(2)+pad | gy < yl(1)-pad | gy > yl(2)+pad;
  gx(bad) = NaN; gy(bad) = NaN;
  if all(isnan(gx)), return; end
  plot(ax, gx, gy, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);
  xlim(ax, xl); ylim(ax, yl);
  drawn = true;
catch
end

end
