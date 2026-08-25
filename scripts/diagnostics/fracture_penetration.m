function OUT = fracture_penetration(opts)
%FRACTURE_PENETRATION How high do the basal crevasses reach into the column?
%
%   REPEAT-PASS PROCESSING, third application. The same repeat-pass
%   machinery this project built for vertical strain rates - per-pair
%   dh(z) from the interferograms, network inversion over every pair,
%   the reference-invariant tide fit - here asked a LOCAL question at
%   named columns instead of a line-wide one, together with the SLC
%   amplitudes the phase chain normally discards.
%
%   THE QUESTION. The tracked bed shows open keels of 6-8 m at the
%   confirmed basal crevasses, but a keel that size cuts the flexural
%   rigidity by only ~7% - far too little for the factor-two softening the
%   local-E* inversion finds seaward. If the crevasses drive the
%   softening, the FRACTURE and its weakened halo must extend well above
%   the open keel the bed pick can see. This diagnostic asks the phase and
%   amplitude data that question directly.
%
%   TWO INDEPENDENT PROBES, both from data already on disk:
%
%   PHASE (the strain-rate machinery, depth-resolved). The vvel pair
%   products carry dh(z) per along-track block - the same arrays the
%   strain-rate inversions use - so the network inversion and the
%   reference-invariant tide fit run PER DEPTH BIN, giving the tidal
%   admittance as a function of depth, adm(z), at every block. Above an
%   open, water-pressurised fracture the column is locally compliant:
%   adm(z) at the crevasse blocks should depart from the flanking control
%   blocks over the depth range the fracture spans, and the shallowest
%   depth where the contrast exceeds its uncertainty bounds the FRACTURE
%   TIP from above.
%
%   AMPLITUDE (the multipass SLCs). Fractured ice scatters: the damaged
%   zone above a keel returns englacial power that intact ice at the same
%   depth does not. The mean power profile at the crevasse columns against
%   the controls localises the scattering zone; its top is a second,
%   independent tip estimate. And the BASAL peak power per pass, regressed
%   on that pass's tide, tests seawater pumping - a water-filled crevasse
%   brightens as the tide loads it.
%
%   TARGETS. The confirmed keel positions, in the driver's x_sea frame
%   (origin at the centre of the first valid 500 m admittance block).
%   The detector itself was removed from the repository by decision on
%   2026-08-23 - it was a partial segmentation - so the positions enter
%   here as constants with their provenance: cross-segment-confirmed
%   incisions from the tracked-bed analysis, final state in git history
%   at 611ed62. Keels: 1.73 km (6.6 m, 4/4 segments), 2.38-2.43 km
%   (7.6/5.9 m, 3-4/4), smaller features 2.64 and 3.19 km. The
%   grounding-zone relief field (0.0-0.6 km) is deliberately NOT a
%   target: slope-break topography is an equally available reading there.
%
%   GL2 is excluded as everywhere in the elasticity chain (its a(x)
%   misfits on the known ref_z anomaly); GL1, GL3, GL4 are analysed.
%
%   opts, all optional:
%     .out_dir    where the png goes
%     .net_dir    all-pairs vvel products
%     .mp_dir     multipass products (for the SLC amplitude slices)
%     .fit_mat    flexure_fit.mat, source of each line's x_sea frame
%     .half_win   half-width of a target/control column window [m]
%                 (default 150 - spans the 30-60 m keels with margin
%                 while staying well inside the 500 m block)
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch \
%       "addpath('<code>/scripts/diagnostics'); fracture_penetration"

if nargin < 1 || isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here)));                       % +vdef

if ~isfield(opts,'out_dir') || isempty(opts.out_dir)
  opts.out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures_flexure';
end
if ~isfield(opts,'net_dir') || isempty(opts.net_dir)
  opts.net_dir = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel_net';
end
if ~isfield(opts,'mp_dir') || isempty(opts.mp_dir)
  opts.mp_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
end
if ~isfield(opts,'fit_mat') || isempty(opts.fit_mat)
  opts.fit_mat = '/kucresis/scratch/hoffmana_sta/vvel/figures_flexure/flexure_fit.mat';
end
if ~isfield(opts,'half_win') || isempty(opts.half_win), opts.half_win = 150; end

LINES  = {'EAGER_2022_GL1','EAGER_2022_GL3','EAGER_2022_GL4'};
KEELS  = [1730 2405 2640 3190];     % m, driver x_sea; 2.38-2.43 as one
KEEL_H = [6.6 7.6 3.6 5.1];         % open keel heights from the tracked bed
% Controls flank the keels on undisturbed floating ice: away from every
% keel, seaward of the grounding-zone relief field, short of the far end
% where coverage thins.
CTRL   = [1250 1980 2900 3600];
MAX_BASELINE = 10;

F = load(opts.fit_mat);
L = F.OUT.lines;

OUT = struct('line',{},'adm_z',{},'amp',{});
for n = 1:numel(LINES)
  pn = LINES{n};
  i  = find(strcmp({L.name}, pn), 1);
  assert(~isempty(i), '%s not in %s', pn, opts.fit_mat);
  sgn = L(i).x_flip_sign; ref = L(i).x_flip_ref;

  %% ---- PHASE: depth-resolved admittance at keel vs control blocks
  f = dir(fullfile(opts.net_dir, [pn '_vvel_*.mat']));
  keep = ~cellfun('isempty', regexp({f.name}, ...
    ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
  f = f(keep);
  assert(~isempty(f), 'no vvel products for %s', pn);

  M = load(fullfile(opts.mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
  Np = numel(M.pass); elev = nan(1,Np); ptime = nan(1,Np);
  for k = 1:Np
    elev(k)  = mean(M.pass(k).elev,'omitnan');
    ptime(k) = mean(M.pass(k).gps_time,'omitnan');
  end

  P = []; W = []; D3 = []; depth = []; blk_xsea = [];
  for q = 1:numel(f)
    tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
      'tokens','once');
    o = load(fullfile(opts.net_dir, f(q).name));
    if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
    if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
    if isempty(depth)
      depth = o.depth_blk(:,1);
      blk_xsea = sgn * (o.Along_track(:) - ref);
      Nz = numel(depth); Nb = numel(blk_xsea);
      D3 = nan(Nz, Nb, 0);
    end
    D3(:,:,end+1) = o.dh_blk; %#ok<AGROW>
    P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
    W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
  end
  assert(size(P,1) >= 10, '%s: only %d usable pairs', pn, size(P,1));

  tday = (ptime - min(ptime))/86400; tide = elev - mean(elev);
  spacing = median(diff(sort(blk_xsea(isfinite(blk_xsea)))));
  fprintf('\n===== %s: %d pairs, %d blocks at %.0f m spacing =====\n', ...
    pn, size(P,1), numel(blk_xsea), spacing);

  % THE FOCUSED LOOK: every block in the keel corridor gets its own
  % depth-resolved admittance - the network inversion runs multirow (depth
  % bins are independent rows), one call per block - so the keels are
  % judged against their own along-track context rather than against a
  % handful of far-away control columns, and the depth onset is structure
  % in a section, not a row of stars in a table.
  WINX = [1200 3600];
  sel_b = find(blk_xsea >= WINX(1) & blk_xsea <= WINX(2));
  [~, ord] = sort(blk_xsea(sel_b)); sel_b = sel_b(ord);
  Nzb = numel(depth); Nsb = numel(sel_b);
  ADM = nan(Nzb, Nsb); SD = nan(Nzb, Nsb);
  for j = 1:Nsb
    rows = squeeze(D3(:, sel_b(j), :));           % Nz x Npair
    N = vdef.invertNetwork(P, rows, struct('n_sigma',3,'weights',W,'n_epoch',Np));
    T = vdef.fitTideAdmittance(N.x, tday, tide);
    ADM(:,j) = T.admittance(:);
    SD(:,j)  = T.admittance_std(:);
  end
  bx = blk_xsea(sel_b);

  % Per keel: its own block against the CORRIDOR MEDIAN of every block
  % that is at least 250 m from every keel. The earlier 300-700 m flanking
  % rule left one or two controls per keel - the keels are only 240-670 m
  % apart, so almost every candidate flank sits inside another keel's
  % exclusion - and a single odd control block then manufactures a whole
  % contrast profile. The median over the ~half-dozen surviving blocks is
  % robust to one of them being weird, and its MAD-based standard error
  % says how well the corridor background is actually known.
  keel_b = arrayfun(@(x) nearest_local(bx, x, 0.55*spacing), KEELS);
  d_any = inf(size(bx));
  for q = 1:numel(KEELS)
    d_any = min(d_any, abs(bx - KEELS(q)));
  end
  cb = find(d_any > 250);
  assert(numel(cb) >= 3, '%s: only %d control blocks in the corridor', pn, numel(cb));
  ctrl_med = median(ADM(:,cb), 2, 'omitnan');
  nctl = sum(isfinite(ADM(:,cb)), 2);
  ctrl_se = 1.4826 * mad(ADM(:,cb), 1, 2) ./ sqrt(max(nctl, 1));
  fprintf('corridor controls: %d blocks (%s km)\n', numel(cb), ...
    strjoin(arrayfun(@(v) sprintf('%.2f',v/1e3), bx(cb), 'uni', 0), ' '));

  A = struct('depth', depth, 'bx', bx, 'ADM', ADM, 'SD', SD, ...
    'keel', nan(Nzb, numel(KEELS)), 'keel_sd', nan(Nzb, numel(KEELS)), ...
    'ctrl_mean', repmat(ctrl_med, 1, numel(KEELS)), ...
    'ctrl_se',   repmat(ctrl_se,  1, numel(KEELS)));
  for j = 1:numel(KEELS)
    if ~isfinite(keel_b(j)), continue; end
    A.keel(:,j)    = ADM(:,keel_b(j));
    A.keel_sd(:,j) = SD(:,keel_b(j));
    fprintf('keel %.2f km -> block %.2f km\n', KEELS(j)/1e3, bx(keel_b(j))/1e3);
  end

  fprintf('depth-resolved contrast, keel block minus LOCAL controls (mm per m tide):\n');
  fprintf('%8s |', 'depth m');
  for j = 1:numel(KEELS), fprintf(' keel %4.2fkm |', KEELS(j)/1e3); end
  fprintf('\n');
  zshow = 40:20:260;
  for z = zshow
    [~, iz] = min(abs(depth - z));
    if ~isfinite(depth(iz)), continue; end
    fprintf('%8.0f |', depth(iz));
    for j = 1:numel(KEELS)
      c  = 1e3*(A.keel(iz,j) - A.ctrl_mean(iz,j));
      cs = 1e3*sqrt(A.keel_sd(iz,j)^2 + A.ctrl_se(iz,j)^2);
      star = ' '; if isfinite(c) && abs(c) > 2*cs, star = '*'; end
      fprintf(' %+6.1f+/-%4.1f%s|', c, cs, star);
    end
    fprintf('\n');
  end
  fprintf('* = |contrast| > 2 sigma. The shallowest starred depth at a keel\n');
  fprintf('bounds that fracture tip from above.\n');

  %% ---- AMPLITUDE: SLC power at keel vs control columns, and vs tide
  amp = amp_slices(pn, opts, sgn, ref, KEELS, CTRL, elev);

  OUT(end+1) = struct('line', pn, 'adm_z', A, 'amp', amp); %#ok<AGROW>
  clear D3 M;
end

%% ---- Pooled across lines: three independent pair sets, one keel each
% The lines are independent networks over near-identical ice, so the
% keel-minus-corridor contrasts pool with inverse-variance weights for a
% sqrt(3)-ish gain - and, more importantly, a keel signature has to be a
% CONSISTENT sign across lines to survive the pooling at all.
fprintf('\n===== POOLED across %d lines: keel minus corridor median =====\n', numel(OUT));
fprintf('%8s |', 'depth m');
for j = 1:numel(KEELS), fprintf(' keel %4.2fkm |', KEELS(j)/1e3); end
fprintf('\n');
depth0 = OUT(1).adm_z.depth;
POOL = nan(numel(depth0), numel(KEELS)); POOLSE = POOL;
for j = 1:numel(KEELS)
  num = zeros(numel(depth0),1); den = zeros(numel(depth0),1);
  for n = 1:numel(OUT)
    A = OUT(n).adm_z;
    c  = A.keel(:,j) - A.ctrl_mean(:,j);
    v  = A.keel_sd(:,j).^2 + A.ctrl_se(:,j).^2;
    ok = isfinite(c) & isfinite(v) & v > 0;
    ci = interp1(A.depth(ok), c(ok)./v(ok), depth0, 'linear', NaN);
    vi = interp1(A.depth(ok), 1./v(ok),    depth0, 'linear', NaN);
    add = isfinite(ci) & isfinite(vi);
    num(add) = num(add) + ci(add);
    den(add) = den(add) + vi(add);
  end
  POOL(den>0, j)   = num(den>0) ./ den(den>0);
  POOLSE(den>0, j) = sqrt(1 ./ den(den>0));
end
for z = 40:20:260
  [~, iz] = min(abs(depth0 - z));
  if ~isfinite(depth0(iz)), continue; end
  fprintf('%8.0f |', depth0(iz));
  for j = 1:numel(KEELS)
    c = 1e3*POOL(iz,j); cs = 1e3*POOLSE(iz,j);
    star = ' '; if isfinite(c) && abs(c) > 2*cs, star = '*'; end
    fprintf(' %+6.1f+/-%4.1f%s|', c, cs, star);
  end
  fprintf('\n');
end
fprintf('* = |pooled contrast| > 2 sigma\n');

%% ---- Figure: the best-confirmed in-band keel (2.38-2.43 km), all lines
draw_summary(OUT, KEELS, KEEL_H, opts);

end

%% ========================================================================
function b = nearest_local(xsea, x, gate)
%NEAREST_LOCAL Index of the block containing x, NaN when none is within
%   just over half a block spacing - so a keel between blocks maps to the
%   one it actually loads, and a keel off the axis maps to nothing.
[d, b] = min(abs(xsea - x));
if d > gate, b = NaN; end
end

%% ========================================================================
function amp = amp_slices(pn, opts, sgn, ref, KEELS, CTRL, elev)
%AMP_SLICES Mean SLC power profiles and basal-power tide regressions.
%   Reads only the column windows it needs from each multipass product
%   (matfile partial reads; the full data cube is 1-1.5 GB), averages
%   power incoherently over the window per pass, and reports:
%     .depth, .keel_dB(z, keel), .ctrl_dB(z)  - relative power profiles
%     .base_tide_slope(keel), _se             - dB per metre of tide at
%                                               the basal peak
%   Depth conversion through the firn column; the basal peak is found per
%   window from the mean profile.
mf = matfile(fullfile(opts.mp_dir, sprintf('%s_multipass03.mat', pn)));
info = whos(mf, 'data');
Nz = info.size(1); Npass = info.size(3);
S = load(fullfile(opts.mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass', 'param_multipass');
mi = S.param_multipass.multipass.baseline_master_idx;
al = S.pass(mi).along_track(:);
% Fast-time axis: the multipass product time base of the main pass
tax = S.pass(mi).time(:);
clear S;

xsea_col = sgn * (al - ref);
P = vdef.firnColumn(vdef.defaultParams());
% depth for twtt BELOW THE SURFACE; surface sits at twtt ~ 0 in these
% products (the project's standing geometry)
dvec = vdef.depthFromTwtt(P, max(tax, 0));

win = opts.half_win;
tide = elev(:) - mean(elev);

grpx = [KEELS CTRL];
nk = numel(KEELS); ng = numel(grpx);
prof = nan(Nz, ng);
base_slope = nan(1, ng); base_se = nan(1, ng);
for g = 1:ng
  cols = find(abs(xsea_col - grpx(g)) <= win);
  if isempty(cols), continue; end
  c0 = min(cols); c1 = max(cols);
  blkd = mf.data(1:Nz, c0:c1, 1:Npass);          % contiguous partial read
  pw = squeeze(mean(abs(blkd).^2, 2, 'omitnan'));  % Nz x Npass
  prof(:,g) = mean(pw, 2, 'omitnan');
  % SELF-NORMALISE each pass by its own shallow englacial power (30-120 m,
  % far above every keel) before the tide regression. Without this, any
  % per-pass bulk gain that happens to correlate with the tide - and GL4's
  % does, at +14 to +17 dB/m IDENTICALLY at every position, which is a
  % calibration signature, not physics - reads as pumping. With it, only
  % power changes CONFINED to the basal return survive.
  nb = dvec > 30 & dvec < 120 & isfinite(prof(:,g));
  pnorm = mean(pw(nb,:), 1, 'omitnan');
  pwn = pw ./ pnorm;
  % Basal peak: strongest return below 150 m, NaN-safe - a plain
  % max(mp .* mask) returns NaN's index the moment the profile has one.
  mp = prof(:,g);
  mp(~isfinite(mp) | dvec <= 150) = -Inf;
  [~, ib] = max(mp);
  band = max(1, ib-3) : min(Nz, ib+3);
  bp = 10*log10(squeeze(mean(pwn(band,:), 1, 'omitnan'))).';
  ok = isfinite(bp) & isfinite(tide);
  if nnz(ok) >= 5
    X = [ones(nnz(ok),1), tide(ok)];             % tide is a column above
    beta = X \ bp(ok);
    r = bp(ok) - X*beta;
    Cv = (sum(r.^2)/(nnz(ok)-2)) * ((X.'*X) \ eye(2));
    base_slope(g) = beta(2); base_se(g) = sqrt(abs(Cv(2,2)));
  end
end

% The pumping test is the DIFFERENCE against the controls' own slope: the
% controls see the same passes, the same normalisation, and the same
% residual bulk effects, so what survives the subtraction is local to the
% keel's basal return.
cw = 1 ./ base_se(nk+1:end).^2;
ctrl_slope = sum(base_slope(nk+1:end) .* cw, 'omitnan') / sum(cw, 'omitnan');
ctrl_slope_se = sqrt(1 / sum(cw, 'omitnan'));

ctrl_dB = 10*log10(mean(prof(:, nk+1:end), 2, 'omitnan'));
amp = struct('depth', dvec, ...
  'keel_dB', 10*log10(prof(:,1:nk)) - ctrl_dB, ...
  'base_tide_slope', base_slope(1:nk) - ctrl_slope, ...
  'base_tide_se', sqrt(base_se(1:nk).^2 + ctrl_slope_se^2), ...
  'ctrl_slope', ctrl_slope, 'ctrl_slope_se', ctrl_slope_se);

fprintf('%s basal power vs tide, keel MINUS control (dB per m; ctrl %+.2f+/-%.2f):', ...
  short(pn), ctrl_slope, ctrl_slope_se);
for g = 1:nk
  fprintf('  %4.2fkm %+5.2f+/-%4.2f', KEELS(g)/1e3, ...
    amp.base_tide_slope(g), amp.base_tide_se(g));
end
fprintf('\n');
end

%% ========================================================================
function s = short(name)
s = strrep(name, 'EAGER_2022_', '');
end

%% ========================================================================
function draw_summary(OUT, KEELS, KEEL_H, opts)
%DRAW_SUMMARY The focused block view: one adm(z, x) section per line over
%   the keel corridor, keel positions ticked, plus the on-keel ladders
%   with their local-control bands. Sections share one diverging scale so
%   the lines read against each other; cells whose fit uncertainty
%   exceeds SD_CAP are blanked rather than drawn as confident colour.
cols = [0.11 0.42 0.69; 0.20 0.60 0.25; 0.75 0.20 0.30];
SD_CAP = 6e-3;                                   % m/m; blank noisier cells
CLIM   = 6;                                      % mm per m of tide
ndiv = 256; half = ndiv/2;
dneg = [0.698 0.094 0.169]; dmid = [0.941 0.937 0.925]; dpos = [0.165 0.471 0.839];
dmap = [interp1([0 1],[dneg; dmid], linspace(0,1,half)); ...
        interp1([0 1],[dmid; dpos], linspace(0,1,ndiv-half))];

hf = figure('Visible','off','Position',[100 100 1240 640],'Color','w');
set(0,'CurrentFigure',hf);

nL = numel(OUT);
for n = 1:nL
  A = OUT(n).adm_z;
  ax = axes('parent',hf,'Position',[0.06+(n-1)*0.245 0.44 0.215 0.50]);
  C = 1e3 * A.ADM;
  C(A.SD > SD_CAP | ~isfinite(A.SD)) = NaN;
  okz = isfinite(A.depth);
  im = imagesc(ax, A.bx/1e3, A.depth(okz), C(okz,:));
  set(im, 'AlphaData', isfinite(C(okz,:)));
  set(ax,'YDir','reverse'); hold(ax,'on');
  for j = 1:numel(KEELS)
    plot(ax, [1 1]*KEELS(j)/1e3, [0 12], 'k-', 'LineWidth', 2);
  end
  colormap(ax, dmap); caxis(ax, [-CLIM CLIM]);
  ylim(ax, [0 280]); xlim(ax, [min(A.bx) max(A.bx)]/1e3);
  grid(ax,'off'); box(ax,'on');
  title(ax, sprintf('(%c)  %s', 'a'+n-1, strrep(OUT(n).line,'EAGER_2022_','')));
  if n == 1, ylabel(ax,'Depth (m)'); else, set(ax,'YTickLabel',[]); end
  xlabel(ax,'Seaward distance (km)');
end
cb = colorbar(ax, 'Position', [0.80 0.44 0.015 0.50]);
set(get(cb,'ylabel'),'string','Admittance (mm per m tide)');

% (d) on-keel ladders against local controls, one panel, all lines
ax4 = axes('parent',hf,'Position',[0.06 0.08 0.40 0.27]);
hold(ax4,'on'); hleg = []; lleg = {};
for n = 1:nL
  A = OUT(n).adm_z;
  for j = 1:size(A.keel,2)
    c  = 1e3*(A.keel(:,j) - A.ctrl_mean(:,j));
    cs = 1e3*sqrt(A.keel_sd(:,j).^2 + A.ctrl_se(:,j).^2);
    ok = isfinite(c) & isfinite(cs) & cs < 1e3*SD_CAP;
    if nnz(ok) < 5, continue; end
    hh = plot(ax4, A.depth(ok), c(ok), '-', 'Color', cols(n,:), 'LineWidth', 1.0);
    sig = ok & abs(c) > 2*cs;
    plot(ax4, A.depth(sig), c(sig), 'o', 'Color', cols(n,:), ...
      'MarkerSize', 4, 'MarkerFaceColor', cols(n,:));
    if j == 1
      hleg(end+1) = hh; lleg{end+1} = strrep(OUT(n).line,'EAGER_2022_',''); %#ok<AGROW>
    end
  end
end
plot(ax4, [0 280], [0 0], 'k-', 'LineWidth', 0.8);
grid(ax4,'on'); box(ax4,'on'); xlim(ax4,[0 280]);
xlabel(ax4,'Depth (m)'); ylabel(ax4,'Keel - local ctrl (mm/m)');
title(ax4,'(d)  On-keel ladders; dots where |contrast| > 2\sigma');
if ~isempty(hleg), legend(ax4, hleg, lleg, 'Location','NorthWest'); end

% (e) the pumping test, keel minus control
ax5 = axes('parent',hf,'Position',[0.55 0.08 0.40 0.27]);
hold(ax5,'on');
for n = 1:nL
  am = OUT(n).amp;
  off = (n - 2) * 0.03;                          % separate the lines
  for k = 1:numel(KEELS)
    if ~isfinite(am.base_tide_se(k)), continue; end
    plot(ax5, KEELS(k)/1e3 + off, am.base_tide_slope(k), 'o', ...
      'Color', cols(n,:), 'MarkerSize', 6, 'MarkerFaceColor', cols(n,:));
    plot(ax5, [1 1]*(KEELS(k)/1e3 + off), ...
      am.base_tide_slope(k) + [-1 1]*am.base_tide_se(k), '-', ...
      'Color', cols(n,:), 'LineWidth', 1.2);
  end
end
plot(ax5, [min(KEELS) max(KEELS)]/1e3 + [-0.2 0.2], [0 0], 'k-', 'LineWidth', 0.8);
grid(ax5,'on'); box(ax5,'on');
xlabel(ax5,'Keel position (km)');
ylabel(ax5,'Basal power vs tide, keel - ctrl (dB/m)');
title(ax5,'(e)  Seawater pumping test');

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end
out_fn = fullfile(opts.out_dir,'EAGER_2022_fracture_penetration.png');
print(hf, out_fn, '-dpng', '-r140');
fprintf('\nWrote %s\n', out_fn);
end
