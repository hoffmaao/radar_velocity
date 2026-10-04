%LINE_SECTIONS Tidal response over the radargram it was measured in, per line.
%
%   FOUR FIGURES, one per repeat-pass line (GL1-GL4), each two panels on a
%   SHARED along-track axis:
%
%     top     mm of column change per metre of tide over the top REF_DEPTH
%             metres - one point per 125 m block with its 1-sigma bar, and
%             the moving-window profile through them, COLOURED on the same
%             diverging scale the map uses and washed toward neutral where
%             the response is under NSIG sigma. A colour here and a colour
%             on the map mean the same number, so a feature can be carried
%             between the two figures by eye
%     bottom  the radargram: multi-pass incoherent mean power in dB against
%             depth, so the internal stratigraphy and the ice base are
%             visible under the measurement made from them
%
%   WHY THE RADARGRAM. It shows the ice the measurement is actually made
%   in - the firn and the englacial stratigraphy of the same 0-REF_DEPTH
%   column the panel above reports - so a feature in the section can be
%   read straight up into the response.
%
%   NOTE THE DEPTH WINDOW IS THE TOP OF THE COLUMN, not the whole ice
%   thickness. That is deliberate: the response is a 0-REF_DEPTH mean, and
%   a section running to the bed put the measured column in the top third
%   of the panel. The consequence is that the ice base near 295 m, and the
%   basal keels an earlier version of this figure was built to show, are
%   NOT in view - DEPTH_MAX brings them back if that is what is wanted.
%
%   THE TWO PANELS ARE THE SAME ICE. The top panel's blocks are averages
%   over the column drawn beneath them, at the same along-track
%   coordinate, so a feature in the section can be read straight up into
%   the measurement. That is the whole point of the pairing, and it is why
%   the x axes are locked rather than merely similar.
%
%   SMOOTHED AT WIN_M, THE SAME AS THE MAP. The response curve is a
%   moving window along the leg - inverse-variance weights times a tricube
%   kernel in along-track distance. The kernel width is the honest
%   resolution; STEP_M only subdivides it for drawing, and neighbouring
%   samples share most of their blocks. The width is deliberately the same
%   as tidal_response_map's so a colour means the same number, smoothed the
%   same way, in both figures.
%
%   The floor on any of this is the PRODUCT: 125 m blocks from
%   CSARP_vvel_net_fine. A kernel narrower than that spacing cannot average
%   anything, and run_vvel_scratch.m's block_size_override is what changes
%   it.
%
%   RADARGRAM CONSTRUCTION. The multipass product is already
%   surface-referenced - pass.surface spans under 3 ns across a whole line,
%   about 0.25 m in ice - so no per-trace flattening is applied and the
%   depth axis comes straight from the fast-time axis through the firn
%   column. Power is averaged INCOHERENTLY over all enabled passes (13-15
%   looks) and then over N_LOOK along-track samples: the passes are
%   coregistered, so this is a legitimate multilook, and it is what makes
%   a basal notch visible against speckle. Residual pass-to-pass
%   misalignment is up to ~2 range bins (~0.55 m), which blurs nothing at
%   the scale this panel is read at.
%
%   NOTE ON THE COLUMN AXES. `data` holds only ENABLED passes while `pass`
%   stays complete, so a data slice k is pass(pass_en_idxs(k)) - the two
%   index spaces diverge as soon as pass_en_mask has a false. `data` also
%   carries a few more columns than the geometry fields (1928 vs 1925 on
%   GL3), so both are truncated to the shorter one.
%
%   Run on the server (needs CSARP_vvel_net_fine):
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../line_sections.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef
addpath(fileparts(mfilename('fullpath')));            % grl_figure
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'diagnostics'));  % pass_tide

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net_fine');      % 125 m blocks
alt_dir = fullfile(root,'CSARP_vvel_net');           % 500 m, per-line fallback
gis_dir = '/kucresis/scratch/hoffmana_sta/vvel/gis';
out_dir = vdef.figureDir();

PASS_NAMES = {'EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
REF_DEPTH = 100; MAX_BASELINE = 10;
WIN_M  = 500;    % moving-window kernel FULL width - the real resolution.
                 % MATCHES tidal_response_map's, so a colour in the section
                 % and a colour on the map mean the same number smoothed the
                 % same way. Changing one without the other silently breaks
                 % that.
STEP_M = 25;     % profile sampling step; rendering only, not resolution
DEPTH_MAX = 400; % m below the surface - the FULL column through the ice
                 % base near 295 m and into the sub-bed returns. The
                 % response above is a 0-REF_DEPTH mean, marked by the
                 % dashed line, so the panel deliberately shows more ice
                 % than the measurement covers: the bed and its keels are
                 % the reason to draw a section at all.
N_LOOK  = 5;     % along-track samples averaged into the radargram
DB_LO = 8; DB_HI = 99.5;    % percentile stretch for the section grayscale,
                 % set for a FULL-column window: over 0-400 m the power
                 % spans a wide range and these percentiles keep the bed
                 % black without crushing the englacial layering. A shallow
                 % window needs wider percentiles, or the thin range it
                 % spans gets mapped across the whole ramp as speckle.
CLIM = 4;        % mm/m colour range for the response curve - the SAME scale
                 % the map uses, so a colour here and a colour there mean
                 % the same number
NSIG = 2;        % sigma for full saturation; below it the curve washes out

PAL.ink = [0 0 0]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.acc = [0.165 0.471 0.839];
% the map's diverging ramp, built identically so the two figures key to
% each other
PAL.div_neg = [0.698 0.094 0.169]; PAL.div_mid = [0.941 0.937 0.925];
PAL.div_pos = [0.165 0.471 0.839];
NDIV = 256;
DMAP = [interp1([0 1],[PAL.div_neg; PAL.div_mid], linspace(0,1,NDIV/2)); ...
        interp1([0 1],[PAL.div_mid; PAL.div_pos], linspace(0,1,NDIV/2))];

Pfirn = vdef.firnColumn(vdef.defaultParams());

for n = 1:numel(PASS_NAMES)
  pn = PASS_NAMES{n};
  fprintf('\n===================== %s =====================\n', pn);

  S = one_line(pn, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE);
  if isempty(S)
    S = one_line(pn, alt_dir, mp_dir, REF_DEPTH, MAX_BASELINE);
    if isempty(S)
      fprintf('*** no usable product in either build - skipped\n');
      continue;
    end
    fprintf('*** absent from the 125 m build, using the 500 m one\n');
  end
  fprintf('%3d blocks at %.0f m, adm %+.2f..%+.2f, median sigma %.2f\n', ...
    numel(S.adm), S.step, min(S.adm), max(S.adm), ...
    median(S.adm_std(isfinite(S.adm_std))));
  min_n = max(2, min(3, floor(WIN_M/max(S.step,eps))));
  [xg, vg, sg, ng] = kernel_profile(S.along, S.adm, S.adm_std, ...
    WIN_M, STEP_M, min_n);
  fprintf('profile %d samples, %d..%d blocks per %.0f m window\n', ...
    nnz(isfinite(vg)), min(ng(ng>0)), max(ng), WIN_M);

  G = radargram(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), ...
    Pfirn, DEPTH_MAX, N_LOOK);
  if isempty(G)
    fprintf('*** radargram unavailable - skipped\n');
    continue;
  end

  draw_line(pn, S, xg, vg, sg, G, REF_DEPTH, WIN_M, DB_LO, DB_HI, ...
    PAL, DMAP, CLIM, NSIG, out_dir, gis_dir);
end

%% ========================================================================
function draw_line(pn, S, xg, vg, sg, G, REF_DEPTH, WIN_M, DB_LO, DB_HI, ...
                   PAL, DMAP, CLIM, NSIG, out_dir, gis_dir)
%DRAW_LINE The two stacked panels for one line.
[h, GRL] = grl_figure(170, 117.5);
set(0,'CurrentFigure',h);
% The panels share a LEFT EDGE and a WIDTH so the two x axes line up
% pixel for pixel; anything less and reading a section feature up into
% the measurement above it is guesswork.
% PW leaves room to the right for the section colour bar AND its tick
% labels: at 0.885 the bar landed at 0.982-0.996 and every label fell off
% the canvas.
PX = 0.075; PW = 0.845;
POS_TOP = [PX 0.615 PW 0.335];
POS_BOT = [PX 0.100 PW 0.455];

xl = [min(G.x) max(G.x)]/1e3;

% ---- bottom: radargram
axb = axes('parent',h,'Position',POS_BOT);
set(0,'CurrentFigure',h);
v = sort(G.db(isfinite(G.db)));
lo = v(max(1,round(DB_LO/100*numel(v))));
hi = v(min(numel(v),round(DB_HI/100*numel(v))));
imagesc(axb, G.x/1e3, G.z, G.db, [lo hi]);
colormap(axb, flipud(gray(256)));
set(axb,'YDir','normal');
axis(axb,'xy'); set(axb,'YDir','reverse');       % depth increases downward
hold(axb,'on');
xlim(axb, xl); ylim(axb, [0 max(G.z)]);
% the depth the response above is measured over
plot(axb, xl, [REF_DEPTH REF_DEPTH], '--', 'Color', PAL.acc, 'LineWidth', 1.4);
set(axb,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off','Layer','top');
xlabel(axb,'Along track (km)','Color',PAL.ink);
ylabel(axb,'Depth (m)','Color',PAL.ink);

% ---- top: the measurement
axt = axes('parent',h,'Position',POS_TOP);
set(0,'CurrentFigure',h); hold(axt,'on');
plot(axt, xl, [0 0], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
ok = isfinite(vg) & isfinite(sg);
if any(ok)
  fill(axt, [xg(ok); flipud(xg(ok))]/1e3, ...
    [vg(ok)-sg(ok); flipud(vg(ok)+sg(ok))], PAL.acc, ...
    'FaceAlpha', 0.18, 'EdgeColor','none');
end
b = isfinite(S.adm) & isfinite(S.adm_std);
ab = S.along(b); vb = S.adm(b); sb = S.adm_std(b);

% Symmetric about zero so the sign of the response reads off the axis,
% and set from the PROFILE rather than the raw block scatter, which on
% the noisier lines runs three times as far.
yr = max([abs(vb); 4]);
yr = min(yr, max(4, 1.6*max(abs(vg(ok)))));
ylim(axt, [-yr yr]); xlim(axt, xl);

ins = abs(vb) <= yr;
for k = 1:numel(ab)
  plot(axt, [ab(k) ab(k)]/1e3, vb(k)+sb(k)*[-1 1], '-', ...
    'Color', PAL.ink_soft, 'LineWidth', 0.7);
end
plot(axt, ab(ins)/1e3, vb(ins), 'o', 'MarkerSize', 4.5, ...
  'MarkerFaceColor', PAL.ink, 'MarkerEdgeColor','w', 'LineWidth', 0.5);
% Blocks off the top or bottom get a hollow chevron ON the boundary. Left
% bare they showed as an error bar running off the panel with no marker,
% which reads as a rendering fault rather than as a real, large value.
ofs = 0.03*yr;
lo = vb < -yr; hi = vb > yr;
if any(lo)
  plot(axt, ab(lo)/1e3, repmat(-yr+ofs, nnz(lo), 1), 'v', 'MarkerSize', 5, ...
    'MarkerFaceColor','none', 'MarkerEdgeColor', PAL.ink, 'LineWidth', 1);
end
if any(hi)
  plot(axt, ab(hi)/1e3, repmat(yr-ofs, nnz(hi), 1), '^', 'MarkerSize', 5, ...
    'MarkerFaceColor','none', 'MarkerEdgeColor', PAL.ink, 'LineWidth', 1);
end
% The response curve, coloured on the map's scale with the map's
% significance rule, so this panel and the map read as one encoding.
% Height already gives the value; the colour adds where it is resolved.
nsig = ribbon_line(axt, xg(ok)/1e3, vg(ok), vg(ok), sg(ok), ...
  DMAP, CLIM, NSIG, PAL.div_mid, 2.6);
fprintf('response curve: %d of %d segments reach %.0f sigma\n', ...
  nsig, max(nnz(ok)-1,0), NSIG);
grid(axt,'on');
set(axt,'GridAlpha',0.15,'XColor',PAL.ink,'YColor',PAL.ink,'Box','off', ...
  'XTickLabel',[]);
ylabel(axt, 'Tidal response (mm/m)','Color',PAL.ink);
if any(lo|hi)
  fprintf('%d of %d blocks on the panel boundary (|adm| > %.1f): %s\n', ...
    nnz(lo|hi), numel(vb), yr, strjoin(arrayfun(@(z) sprintf('%+.1f', z), ...
    vb(lo|hi).', 'UniformOutput', false), ', '));
end

% grounding-line crossing, on BOTH panels
xgl = gl_crossing(S, gis_dir);
if isfinite(xgl)
  plot(axt, [xgl xgl]/1e3, ylim(axt), ':', 'Color', PAL.ink, 'LineWidth', 1.2);
  plot(axb, [xgl xgl]/1e3, ylim(axb), ':', 'Color', PAL.ink, 'LineWidth', 1.2);
  fprintf('MEaSUREs grounding line crosses at %.2f km along\n', xgl/1e3);
end

% no panel titles: lettering and captions go in the slide or manuscript,
% same convention as the rest of the figures here
cb = colorbar(axb,'eastoutside');
set(axb,'Position',POS_BOT);
set(cb,'Position',[PX+PW+0.012 POS_BOT(2) 0.014 POS_BOT(4)], ...
  'XColor',PAL.ink,'YColor',PAL.ink);
set(cb.Label,'String','Power (dB)','Color',PAL.ink);

nm = strrep(pn,'EAGER_2022_','');
out_fn = fullfile(out_dir, sprintf('EAGER_2022_%s_section.png', nm));
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi));
close(h);
fprintf('Wrote %s (%.0f m kernel)\n', out_fn, WIN_M);
end

%% ========================================================================
function G = radargram(fn, Pfirn, depth_max, n_look)
%RADARGRAM Multi-pass incoherent mean power, in dB, against depth.
%   Reads only the fast-time window that reaches depth_max, through
%   matfile where the file allows it, because the full array is ~850 MB
%   per product and only about a third of it is above the ice base.
G = [];
if ~exist(fn,'file')
  fprintf('multipass product missing: %s\n', fn); return;
end
try
  L = load(fn, 'pass', 'param_multipass');
  t = L.pass(1).time(:);
  at = L.pass(1).along_track(:);

  % fast-time window: surface sits at twtt ~0 in these products
  zt = vdef.depthFromTwtt(Pfirn, max(t, 0));
  b0 = find(t >= 0, 1, 'first');
  b1 = find(zt <= depth_max, 1, 'last');
  if isempty(b0) || isempty(b1) || b1 <= b0
    fprintf('no usable fast-time window\n'); return;
  end

  D = [];
  try
    m = matfile(fn);
    D = m.data(b0:b1, :, :);
  catch
    Q = load(fn, 'data');
    D = Q.data(b0:b1, :, :);
    clear Q;
  end

  % Columns: `data` can carry a few more than the geometry fields
  nc = min(size(D,2), numel(at));
  D = D(:, 1:nc, :); at = at(1:nc);

  P = mean(abs(double(D)).^2, 3, 'omitnan');     % incoherent, over passes
  clear D;
  if n_look > 1
    P = movmean(P, n_look, 2, 'omitnan');        % and along track
  end
  db = 10*log10(P);
  db = db - max(db(isfinite(db)));               % relative to the peak

  G = struct('x', at, 'z', zt(b0:b1), 'db', db, ...
    'npass', size(L.pass,2), 'bins', [b0 b1]);
  fprintf('radargram %d bins (0-%.0f m) x %d cols, %d passes, %d-look\n', ...
    b1-b0+1, max(G.z), nc, numel(L.pass), n_look);
catch ME
  fprintf('radargram failed (%s)\n', ME.message);
  G = [];
end
end

%% ========================================================================
function xgl = gl_crossing(S, gis_dir)
%GL_CROSSING Along-track position where the leg crosses the MEaSUREs line.
%   The track actually crosses the polyline, so the crossing is where the
%   distance to it is minimised; NaN (and no annotation) if the nearest
%   approach is too far to call it a crossing.
%
%   The polyline is DENSIFIED first. MEaSUREs vertices are 0.5-2 km apart
%   in this window, so nearest-VERTEX distance stays in the hundreds of
%   metres even where the track crosses the line dead on - measuring to
%   vertices found no crossing on any of the four legs.
xgl = NaN;
fn = fullfile(gis_dir,'GroundingLine_Antarctica_v02.shp');
if ~exist(fn,'file'), return; end
try
  ps = projcrs(3031);
  [tx, ty] = projfwd(ps, S.track_lat, S.track_lon);
  Q = shaperead(fn);
  gx = []; gy = [];
  for k = 1:numel(Q)
    gx = [gx; Q(k).X(:)]; gy = [gy; Q(k).Y(:)]; %#ok<AGROW>
  end
  near = gx > min(tx)-5e3 & gx < max(tx)+5e3 & ...
         gy > min(ty)-5e3 & gy < max(ty)+5e3 & isfinite(gx);
  gx = gx(near); gy = gy(near);
  if numel(gx) < 2, return; end
  % densify to ~20 m so a crossing between vertices is actually found
  dx = []; dy = [];
  for k = 1:numel(gx)-1
    L = hypot(gx(k+1)-gx(k), gy(k+1)-gy(k));
    if ~isfinite(L) || L == 0 || L > 5e3, continue; end
    t = linspace(0, 1, max(2, ceil(L/20))).';
    dx = [dx; gx(k) + t*(gx(k+1)-gx(k))]; %#ok<AGROW>
    dy = [dy; gy(k) + t*(gy(k+1)-gy(k))]; %#ok<AGROW>
  end
  if isempty(dx), return; end
  d = inf(numel(tx),1);
  for k = 1:numel(tx)
    d(k) = min(hypot(dx - tx(k), dy - ty(k)));
  end
  [dm, k] = min(d);
  if dm < 150
    xgl = S.track_at(k);
  else
    fprintf('nearest approach to the grounding line is %.0f m - not marked\n', dm);
  end
catch
end
end

%% ========================================================================
function S = one_line(pn, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE)
%ONE_LINE Per-block tidal admittance and the leg's full-resolution track.
S = [];
f = dir(fullfile(net_dir, [pn '_vvel_*.mat']));
if isempty(f), return; end
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f), return; end
L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end

P = []; D = []; W = []; along = []; Nblk = 0; imain = 1;
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  o = load(fullfile(net_dir, f(q).name));
  if ~vdef.pairAligned(o), continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  if Nblk == 0
    Nblk = numel(o.S1); along = o.Along_track(:);
    if isfield(o,'baseline_main_idx') && ~isempty(o.baseline_main_idx)
      imain = o.baseline_main_idx;
    end
  end
  sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
  end
  if all(~isfinite(sv)), continue; end
  if isempty(D), D = sv; else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < 10, return; end

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));
tday = (ptime - min(ptime))/86400;
% The tide is CATS2008 at each pass mid-time, with the pass gate, from the
% same helper the flexure inversion uses (scripts/diagnostics/pass_tide.m)
% - NOT the line-mean GPS height, which carries a per-pass height error
% that attenuated this admittance by 15-25% and let a 1.1 m bad pass
% through (see vdef.surfaceAdmittance).
tide = pass_tide(pn, mp_dir);
A = vdef.fitTideAdmittance(N.x, tday, tide);

if imain < 1 || imain > Np, imain = 1; end
S = struct('name',pn, 'along',along, ...
  'adm',     1e3*REF_DEPTH*A.admittance(:), ...
  'adm_std', 1e3*REF_DEPTH*A.admittance_std(:), ...
  'step',    median(diff(along)), ...
  'track_at',  L.pass(imain).along_track(:), ...
  'track_lat', L.pass(imain).lat(:), ...
  'track_lon', L.pass(imain).lon(:));
end



%% ========================================================================
function [xg, vg, sg, ng] = kernel_profile(along, val, sigma, W, step, min_n)
%KERNEL_PROFILE Locally weighted inverse-variance mean along track.
%   Blocks inside the window are weighted by 1/sigma^2 TIMES a tricube
%   kernel in along-track distance, so they fade in and out rather than
%   stepping in - a boxcar over ~40 blocks per leg puts the window's own
%   edges into the profile as if they were structure.
%
%   W is the FULL kernel width and is the honest along-track resolution.
%   `step` only sets how finely the result is sampled for drawing;
%   neighbouring samples share most of their blocks and are strongly
%   correlated, so `step` must never be quoted as a resolution.
%
%   sg is the standard error of the WEIGHTED MEAN, sqrt(sum(w^2 s^2))/sum(w).
%   It assumes the blocks are independent, which is the same assumption
%   the disjoint-bin stack it replaces was making.
ok = isfinite(along) & isfinite(val) & isfinite(sigma) & sigma > 0;
along = along(ok); val = val(ok); sigma = sigma(ok);
if isempty(along)
  xg = []; vg = []; sg = []; ng = []; return;
end
xg = (min(along) : step : max(along)).';
vg = nan(size(xg)); sg = nan(size(xg)); ng = zeros(size(xg));
for k = 1:numel(xg)
  u  = abs(along - xg(k)) / (W/2);
  in = u < 1;
  if nnz(in) < min_n, continue; end
  kern = (1 - u(in).^3).^3;                 % tricube
  w  = kern ./ sigma(in).^2;
  sw = sum(w);
  if sw <= 0, continue; end
  vg(k) = sum(w .* val(in)) / sw;
  sg(k) = sqrt(sum((w.^2) .* (sigma(in).^2))) / sw;
  ng(k) = nnz(in);
end
end
