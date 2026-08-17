%ERROR_BUDGET What limits this measurement, and what would move it.
%
%   The companion to tidal_evidence.m. That figure shows the tidal result
%   is not detectable; this one shows WHY, by putting every error term
%   that has been measured on one axis against the signal.
%
%   (a) THE BUDGET. Each bar is a term that has actually been measured,
%       not modelled. Two of them are upper bounds that turned out small,
%       which is the useful part: they are eliminations. The signal is
%       drawn as a line so the gap is read directly.
%
%   (b) THE NETWORK'S OWN VERDICT. Closure residual per pair from the
%       full 78-pair inversion, with the pairs it rejected marked. The
%       network identifies its bad data without being told what to look
%       for, and the offenders concentrate on one pass.
%
%   (c) PER-PASS RESIDUAL. The same information collapsed onto passes,
%       which is where it is actionable - one pass dominates.
%
%   (d) WHAT MORE EPOCHS CANNOT FIX. Formal precision falls as
%       1/sqrt(N_passes). At 13 passes the single-line formal precision
%       (~1.5 mm) is comparable to - marginally above - the 1.24 mm ApRES
%       rate-method signal, and pooling the four calibrated lines
%       (0.55 mm, line_means.m) puts the measurement clearly below it. The
%       between-build disagreement is not a precision term and does not
%       move with N at all. Panel (d) draws both, and the gap between them
%       is the honest statement of what is wrong: the formal error badly
%       understates the true one, and the term that dominates has not yet
%       been diagnosed. An earlier reading of this budget called more
%       epochs "the only lever left"; that was wrong, and this panel is
%       what corrects it - more epochs cannot cross the between-build
%       reproducibility floor.
%
%   Requires the network product (pairs='all'):
%     matlab -batch "only_pass_names={'EAGER_2022_GL3'}; pairing_override='all'; out_suffix='_net'; run('.../run_vvel_scratch.m')"
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../error_budget.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net');
out_dir = '/kucresis/scratch/hoffmana_sta/vvel/figures';
PASS_NAME = 'EAGER_2022_GL3';
REF_DEPTH = 100; MAX_BASELINE = 10;

% Measured elsewhere in the project, cited so the bars are traceable
SIG_APRES   = 1.24;   % mm per m of tide over the top 100 m (magnitude of
                      % -1.24 +/- 0.04), ApRES GA04 rate method,
                      % apres_rate_check.py
FLOOR_BUILD = 8.28;   % mm, EAGER_2022 vs GL1 same-leg disagreement
% DECOMPOSED 2026-08-17 (diagnostics/between_build.m, master_sensitivity.m):
% the between-build disagreement quoted against what the two builds' own
% sigmas predict is 1.2-1.5x for calibrated pairs and 2.6x for EAGER_2022
% vs GL1. In quadrature-excess terms per block that is a median 3.6 mm
% among GL1-GL4 and 7.7 mm for the EAGER_2022 comparison. The master-pass
% choice is CLEARED: same input, same passes, master 11 vs 6 differs by
% 0.7x expected - i.e. measured zero excess.
XS_CALIB = 3.6;       % mm/block, median quadrature excess, calibrated pairs
XS_EAGER = 7.7;       % mm/block, EAGER_2022 vs GL1 quadrature excess
N_ALLLEGS   = 50;     % traverses an all-legs combine_passes rebuild would give

PAL.cat = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; ...
           0.929 0.631 0.000; 0.910 0.482 0.643];
PAL.ink = [0.20 0.20 0.20]; PAL.ink_soft = [0.45 0.45 0.45];
PAL.band = [0.90 0.90 0.88]; PAL.expect = [0.35 0.35 0.35];

%% Rebuild the network result
L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', PASS_NAME)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k) = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
clear L;

f = dir(fullfile(net_dir, [PASS_NAME '_vvel_*.mat']));
P = []; D = []; W = [];
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',PASS_NAME) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  if isempty(tok), continue; end
  o = load(fullfile(net_dir, f(q).name));
  if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  Nblk = numel(o.S1); sv = nan(Nblk,1);
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
N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));

tday = (ptime - min(ptime))/86400;
tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);
mm = @(v) 1e3*REF_DEPTH*v;

% terms, all in mm
closure_mm = mm(median(N.rms,'omitnan'));
fit_mm     = median(mm(A.admittance_std(isfinite(A.admittance_std))));
% per-pass systematic: how much of the residual repeats block to block
res_ep = nan(size(N.x));
for b = 1:size(N.x,1)
  y = N.x(b,:); ok = isfinite(y) & isfinite(tday) & isfinite(tide);
  if nnz(ok) < 6, continue; end
  X = [ones(nnz(ok),1), tday(ok).'-mean(tday(ok)), tide(ok).'-mean(tide(ok))];
  res_ep(b,ok) = (y(ok).' - X*(X\y(ok).')).';
end
rho = [];
gd = find(sum(isfinite(res_ep),2) >= 6);
for a1 = 1:numel(gd), for a2 = a1+1:numel(gd)
  u = res_ep(gd(a1),:); v = res_ep(gd(a2),:); ok = isfinite(u)&isfinite(v);
  if nnz(ok)>=6 && std(u(ok))>0 && std(v(ok))>0
    m = corrcoef(u(ok),v(ok)); rho(end+1)=m(1,2); %#ok<AGROW>
  end
end, end
rho_med = median(rho);
tot_res = mm(median(std(res_ep,0,2,'omitnan'),'omitnan'));
perpass_mm = tot_res * sqrt(max(rho_med,0));    % the part that repeats
model_mm   = 0.05 * tot_res;                    % phase-lag test: 5% reduction
unexpl_mm  = sqrt(max(tot_res^2 - closure_mm^2 - perpass_mm^2 - model_mm^2, 0));

fprintf('\n=== error budget, mm of column change (top %d m) ===\n', REF_DEPTH);
fprintf('  per-pair noise (closure)      %6.2f\n', closure_mm);
fprintf('  per-pass systematic           %6.2f   (block-to-block rho = %.2f)\n', perpass_mm, rho_med);
fprintf('  model form (tidal phase lag)  %6.2f   (5%% of scatter)\n', model_mm);
fprintf('  unexplained, per block/epoch  %6.2f\n', unexpl_mm);
fprintf('  --- within-product total      %6.2f\n', tot_res);
fprintf('  between-build (same leg)      %6.2f\n', FLOOR_BUILD);
fprintf('  ApRES signal (rate method)    %6.2f\n', SIG_APRES);

%% Figure
h = figure('Visible','off','Position',[100 100 1120 880],'Color','w');
set(0,'CurrentFigure',h);
axst = {'GridAlpha',0.15,'XColor',PAL.ink_soft,'YColor',PAL.ink_soft,'Box','off'};

% (a) budget
ax1 = axes('parent',h,'Position',[0.20 0.585 0.30 0.345]);
set(0,'CurrentFigure',h); hold(ax1,'on');
lbl = {'per-pair noise','per-pass systematic','model (phase lag)', ...
       'master choice (cleared)','unexplained','WITHIN-PRODUCT', ...
       'build excess, calibrated','build excess, EAGER\_2022'};
val = [closure_mm perpass_mm model_mm 0 unexpl_mm tot_res XS_CALIB XS_EAGER];
col = [PAL.cat(3,:); PAL.cat(3,:); PAL.cat(3,:); PAL.cat(3,:); ...
       PAL.cat(2,:); PAL.ink_soft; PAL.cat(2,:); PAL.cat(2,:)];
for k = 1:numel(val)
  barh(ax1, k, val(k), 0.62, 'FaceColor', col(k,:), 'EdgeColor','none');
end
plot(ax1, [SIG_APRES SIG_APRES], [0.4 numel(val)+0.6], '--', ...
  'Color', PAL.ink, 'LineWidth', 2);
text(ax1, SIG_APRES, numel(val)+0.55, ' ApRES (rate method)', 'Color', PAL.ink, ...
  'FontSize', 9, 'VerticalAlignment','top');
set(ax1,'YTick',1:numel(val),'YTickLabel',lbl,'YDir','reverse');
grid(ax1,'on'); set(ax1, axst{:}); ylim(ax1,[0.4 numel(val)+0.6]);
xlabel(ax1,'mm of column change, top 100 m','Color',PAL.ink);
title(ax1,'(a) Everything that has been measured (green = eliminated)','Color',PAL.ink);

% (b) closure per pair
ax2 = axes('parent',h,'Position',[0.60 0.585 0.355 0.345]);
set(0,'CurrentFigure',h); hold(ax2,'on');
% only rows the inversion actually solved can reject a pair; unsolved rows
% leave N.used all-false without meaning rejection
solved = any(isfinite(N.x),2);
rej = sum(~N.used(solved,:) & isfinite(D(solved,:)),1);
cr  = mm(sqrt(mean(N.resid.^2,1,'omitnan')));
kept = rej == 0;
plot(ax2, find(kept), cr(kept), 'o', 'MarkerSize', 6, ...
  'MarkerFaceColor', PAL.cat(1,:), 'MarkerEdgeColor','w','LineWidth',0.8);
plot(ax2, find(~kept), cr(~kept), 'o', 'MarkerSize', 9, ...
  'MarkerFaceColor', PAL.cat(2,:), 'MarkerEdgeColor', PAL.ink, 'LineWidth',1.2);
grid(ax2,'on'); set(ax2, axst{:});
xlabel(ax2,'pair index','Color',PAL.ink);
ylabel(ax2,'closure residual (mm)','Color',PAL.ink);
title(ax2, sprintf('(b) %d of %d pairs rejected in at least one block', ...
  nnz(~kept), numel(kept)), 'Color', PAL.ink);

% (c) per-pass residual
ax3 = axes('parent',h,'Position',[0.20 0.09 0.30 0.345]);
set(0,'CurrentFigure',h); hold(ax3,'on');
rms_ep = mm(sqrt(mean(res_ep.^2,1,'omitnan')));
for k = 1:numel(rms_ep)
  if ~isfinite(rms_ep(k)), continue; end
  cc = PAL.cat(1,:); if rms_ep(k) > 2*median(rms_ep,'omitnan'), cc = PAL.cat(2,:); end
  bar(ax3, k, rms_ep(k), 0.65, 'FaceColor', cc, 'EdgeColor','none');
end
grid(ax3,'on'); set(ax3, axst{:});
xlabel(ax3,'pass','Color',PAL.ink);
ylabel(ax3,'post-fit residual (mm)','Color',PAL.ink);
title(ax3,'(c) One pass dominates','Color',PAL.ink);

% (d) what would fix it
ax4 = axes('parent',h,'Position',[0.60 0.09 0.355 0.345]);
set(0,'CurrentFigure',h); hold(ax4,'on');
nep = 5:2:60;
proj = fit_mm * sqrt(size(N.x,2) ./ nep);       % formal precision ~ 1/sqrt(N)
plot(ax4, nep, proj, '-', 'Color', PAL.cat(1,:), 'LineWidth', 2.2);
% the between-build disagreement is not a precision term: it does not fall
% with N, so it is drawn flat. This is the whole point of the panel.
plot(ax4, [min(nep) max(nep)], [FLOOR_BUILD FLOOR_BUILD], '-', ...
  'Color', PAL.cat(2,:), 'LineWidth', 2.6);
plot(ax4, [min(nep) max(nep)], [SIG_APRES SIG_APRES], '--', ...
  'Color', PAL.ink, 'LineWidth', 2);
plot(ax4, size(N.x,2), fit_mm, 'o', 'MarkerSize', 11, ...
  'MarkerFaceColor', PAL.cat(1,:), 'MarkerEdgeColor','w','LineWidth',1.2);
pr50 = fit_mm*sqrt(size(N.x,2)/N_ALLLEGS);
plot(ax4, N_ALLLEGS, pr50, 'p', 'MarkerSize', 16, ...
  'MarkerFaceColor', PAL.cat(4,:), 'MarkerEdgeColor', PAL.ink, 'LineWidth',1.2);
text(ax4, max(nep), FLOOR_BUILD, 'between-build disagreement ', ...
  'Color', PAL.cat(2,:), 'FontSize', 9, 'HorizontalAlignment','right', ...
  'VerticalAlignment','bottom');
text(ax4, max(nep), SIG_APRES, 'ApRES (rate method) ', 'Color', PAL.ink, ...
  'FontSize', 9, 'HorizontalAlignment','right','VerticalAlignment','bottom');
text(ax4, size(N.x,2), fit_mm, sprintf('  formal, now (%d passes)', size(N.x,2)), ...
  'Color', PAL.ink, 'FontSize', 9, 'VerticalAlignment','top');
text(ax4, N_ALLLEGS, pr50, sprintf('all legs (~%d)  ', N_ALLLEGS), ...
  'Color', PAL.ink, 'FontSize', 9, 'VerticalAlignment','top', ...
  'HorizontalAlignment','right');
grid(ax4,'on'); set(ax4, axst{:}); xlim(ax4,[min(nep) max(nep)]);
ylim(ax4,[0 1.15*FLOOR_BUILD]);
xlabel(ax4,'number of passes','Color',PAL.ink);
ylabel(ax4,'1-sigma (mm)','Color',PAL.ink);
title(ax4,'(d) More epochs cannot cross the reproducibility floor','Color',PAL.ink);

print(h, fullfile(out_dir,'EAGER_2022_error_budget.png'), '-dpng','-r120');
close(h);
fprintf('\nWrote %s\n', fullfile(out_dir,'EAGER_2022_error_budget.png'));
