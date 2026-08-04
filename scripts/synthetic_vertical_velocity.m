%SYNTHETIC_VERTICAL_VELOCITY End-to-end validation of the +vdef chain.
%
%   Builds a synthetic repeat-pass SLC pair from a KNOWN vertical strain
%   rate profile, pushes it through the full processing chain
%
%     vdef.multilook -> vdef.differentialRange -> vdef.blockAverage
%       -> vdef.verticalDisplacement -> vdef.invertStrainRate
%
%   and checks that the inversion recovers the truth. The synthetic is
%   scaled to the McMurdo repeat-pass geometry the project targets:
%   accum3 at 750 MHz, 300 MHz sampling, ~550 m of coherent record,
%   and a repeat interval of a few days.
%
%   The point of the exercise is that the per-pixel signal here is far
%   BELOW the per-pixel phase noise - a few days of firn compaction moves
%   the deepest reflector by about a millimetre, i.e. a few hundredths of a
%   radian - so the test also demonstrates that the precision comes from
%   coherence-weighted along-track averaging, not from any single look.
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui synthetic_vertical_velocity.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef

rand('seed', 7); randn('seed', 7);   %#ok<RAND> % Octave-compatible seeding

C = vdef.constants();

%% 1. Geometry and truth
% =====================================================================
fc      = 750e6;             % accum3 centre frequency [Hz]
fs      = 300e6;             % sampling rate [Hz]
Nt      = 2000;              % ~6.7 us record
Nx      = 4000;              % range lines
Surface = 0.20e-6;           % surface twtt [s]
delta_t = 3 / 365.25;        % repeat interval [yr]

Time = (0:Nt-1).' / fs;

par = vdef.defaultParams();
par.bco_depth = 60;

% Truth: vertical strain rate linear in depth over the fitted range.
% Negative = vertical compression (the column is shortening).
H_norm = 500;
S1_true = -1.2e-3;           % depth-averaged vertical strain rate [1/yr]
S2_true =  0.6e-3;           % top-to-bottom gradient [1/yr]
eps_zz_true = @(d) S1_true + S2_true*(2*d/H_norm - 1);

F = vdef.forwardDisplacement(eps_zz_true, Time, Surface, par, delta_t);
dtau_true = F.dtau(:,1);

fprintf('Truth: S1 = %+.3e /yr, S2 = %+.3e /yr\n', S1_true, S2_true);
fprintf('Peak |dtau| over the record: %.1f ps (%.3f rad at %.0f MHz)\n', ...
  max(abs(dtau_true))*1e12, 2*pi*fc*max(abs(dtau_true)), fc/1e6);

%% 2. Synthetic SLC pair
% =====================================================================
% ref is a fully developed speckle field with a depth-decaying envelope.
% sec is the same field delayed by dtau, plus a decorrelation component
% whose weight grows with depth (coherence loss with range).
depth_col = F.depth(:,1);
amp = 10.^(-depth_col/1200);            % power roll-off with depth
amp(~isfinite(amp)) = 0;

gamma = 0.90 * exp(-max(depth_col,0)/700);   % true coherence vs depth
gamma(~isfinite(gamma)) = 0;

s_ref = bsxfun(@times, amp, (randn(Nt,Nx) + 1i*randn(Nt,Nx))/sqrt(2));
indep = bsxfun(@times, amp, (randn(Nt,Nx) + 1i*randn(Nt,Nx))/sqrt(2));

% Matched-filter convention: a delay dtau multiplies the signal by
% exp(-1i*2*pi*fc*dtau), so sec.*conj(ref) has phase -2*pi*fc*dtau and the
% correct recovery sign is phase_sign = -1.
carrier = exp(-1i*2*pi*fc*dtau_true);
carrier(~isfinite(carrier)) = 1;

s_sec = bsxfun(@times, sqrt(gamma) .* carrier, s_ref) ...
      + bsxfun(@times, sqrt(max(1-gamma,0)), indep);

%% 3. Processing chain
% =====================================================================
opts = [];
opts.mlook_window        = [5 15];
opts.phase_sign          = -1;
opts.ref_twtt_offset     = 50e-9;
opts.coherence_threshold = 0.3;
opts.max_gap_bins        = 20;
opts.min_coverage        = 0.3;
opts.block_size          = 2000;
opts.order               = 2;
opts.fit_top_depth       = 20;
opts.fit_bot_depth       = 500;
opts.norm_depth          = H_norm;
opts.reg                 = 0;
opts.bins_per_look       = opts.mlook_window(1);
opts.cols_per_look       = opts.mlook_window(2);
opts.min_samples         = 50;
opts.delta_t             = delta_t;
opts.densification_rate  = 0;

[igram, coh] = vdef.multilook(s_ref, s_sec, opts.mlook_window);

map = [];
map.Time      = Time;
map.Surface   = repmat(Surface, 1, Nx);
map.fc        = fc;
map.phase     = angle(igram);
map.coherence = coh;
map.phase_is_unwrapped = false;

[dtau, info] = vdef.differentialRange(map, opts);
blk = vdef.blockAverage(dtau, map, info, opts);
V   = vdef.verticalDisplacement(blk, map, par, opts);
S   = vdef.invertStrainRate(V, blk, opts);

%% 4. Report
% =====================================================================
fprintf('\nBlocks inverted: %d\n', nnz(isfinite(S.S1)));
fprintf('%-6s %12s %12s %12s %10s %10s\n', 'block', 'S1 [1/yr]', 'S2 [1/yr]', 'rms [m/yr]', 'p_quad', 'n_eff');
for b = 1:numel(S.S1)
  fprintf('%-6d %12.4e %12.4e %12.4e %10.3g %10.0f\n', ...
    b, S.S1(b), S.S2(b), S.rms(b), S.p_quad(b), S.n_eff(b));
end

S1_err = S.S1 - S1_true;
S2_err = S.S2 - S2_true;
fprintf('\nS1 error: %+.3e /yr (%.1f%% of truth)\n', ...
  mean(S1_err,'omitnan'), 100*mean(S1_err,'omitnan')/abs(S1_true));
fprintf('S2 error: %+.3e /yr (%.1f%% of truth)\n', ...
  mean(S2_err,'omitnan'), 100*mean(S2_err,'omitnan')/abs(S2_true));

% Set just outside the errors this seeded synthetic actually produces (S1
% 0.3%, S2 1.0%) rather than at a token 20%, so that a change to the fit
% weighting trips this instead of passing silently. The seed is fixed and
% the numbers reproduce bit-for-bit, so the margin only has to cover
% platform arithmetic, not run-to-run scatter.
tol_S1 = 0.02;
tol_S2 = 0.03;   % the gradient is the weaker-constrained coefficient
assert(abs(mean(S1_err,'omitnan')) < tol_S1*abs(S1_true), ...
  'S1 not recovered within %.0f%%', 100*tol_S1);
assert(abs(mean(S2_err,'omitnan')) < tol_S2*abs(S2_true), ...
  'S2 not recovered within %.0f%%', 100*tol_S2);
fprintf('\nPASS: strain rate recovered within %.0f%%/%.0f%% of truth.\n', ...
  100*tol_S1, 100*tol_S2);

%% 5. Degenerate per-sample sigma
% =====================================================================
% A bin whose within-block scatter collapses to zero, or is missing
% altogether, carries no weighting information. It has to be excluded, not
% given the largest weight in the block: clamping the sigma instead of
% masking it hands that one bin a weight of 1/eps^2 and the normalisation
% then crushes every legitimate sample to ~1e-31, so the fit is driven
% entirely by the degenerate bin.
Vd = V;
fit_rows = find(isfinite(V.depth(:,1)) & isfinite(V.v(:,1)) & ...
  V.depth(:,1) >= opts.fit_top_depth & V.depth(:,1) <= opts.fit_bot_depth);
assert(numel(fit_rows) >= 2, 'no fitted samples to degrade');
Vd.v_scatter(fit_rows(1),   1) = 0;
Vd.v_scatter(fit_rows(end), 1) = NaN;
Sd = vdef.invertStrainRate(Vd, blk, opts);

fprintf('\nBlock 1 with one zero and one NaN sigma:\n');
fprintf('  S1 %+.4e -> %+.4e /yr (truth %+.4e)\n', S.S1(1), Sd.S1(1), S1_true);
fprintf('  S2 %+.4e -> %+.4e /yr (truth %+.4e)\n', S.S2(1), Sd.S2(1), S2_true);

assert(isfinite(Sd.S1(1)) && isfinite(Sd.S2(1)), ...
  'degenerate sigma left block 1 uninverted');
% Block 1 on its own is noisier than the line mean asserted above - its S2
% sits at 3.4% - so this check carries its own, looser tolerance.
tol_blk = 0.07;
assert(abs(Sd.S1(1) - S1_true) < tol_blk*abs(S1_true) && ...
       abs(Sd.S2(1) - S2_true) < tol_blk*abs(S2_true), ...
  'degenerate sigma broke the recovery of block 1');
% Dropping two samples out of a few hundred should barely move the fit
assert(abs(Sd.S1(1) - S.S1(1)) < 0.02*abs(S1_true) && ...
       abs(Sd.S2(1) - S.S2(1)) < 0.02*abs(S2_true), ...
  'degenerate sigma perturbed block 1 far more than dropping two samples should');
fprintf('PASS: degenerate sigma samples are masked, not weighted up.\n');

%% 6. Figure
% =====================================================================
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Three panels side by side do not fit in the default 560 pt figure width:
% gnuplot keeps the tick count and lets the labels run together, so
% '-80-60-40-20' reads as one number. Widening the figure is the portable
% fix - xticks/xtickformat are not available in both MATLAB and Octave.
h = figure('Visible','off','Position',[100 100 1000 420]);

subplot(1,3,1);
plot(1e12*blk.dtau, V.depth(:,1), '-'); hold on;
plot(1e12*dtau_true, depth_col, 'k--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); ylim([0 550]); grid on;
xlabel('\Delta\tau (ps)'); ylabel('Depth (m)');
title('Differential traveltime');
% NorthWest, not SouthWest: dtau grows downward from zero at the surface,
% so the bottom-left corner is where the deep samples are and the legend
% would sit on top of them.
legend('blocks','truth','Location','NorthWest');

subplot(1,3,2);
plot(V.v, V.depth(:,1), '-'); hold on;
plot(S.v_fit, S.depth_grid, 'k-', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); ylim([0 550]); grid on;
xlabel('Relative w (m/yr)'); title('Vertical velocity');

subplot(1,3,3);
% Strain rates here are a few times 1e-3/yr, so plot in units of 1e-3/yr:
% the raw values give tick labels like '-0.0018' that collide into an
% unreadable run under the gnuplot backend. Scaling the data beats
% xticks/xtickformat, which are not portable between MATLAB and Octave.
eps_scale = 1e3;
dd = linspace(0,H_norm,100);
plot(eps_scale*S.eps_zz, S.depth_grid, '-', 'LineWidth', 1.5); hold on;
plot(eps_scale*eps_zz_true(dd), dd, 'k--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); ylim([0 550]); grid on;
% Round the x limits out to whole 1e-3/yr so gnuplot lays down a few
% widely spaced ticks rather than one per 0.2 across a narrow panel.
eps_all = eps_scale*[S.eps_zz(:); eps_zz_true(dd(:))];   % min/max skip NaN
xlim([floor(min(eps_all)) ceil(max(eps_all))]);
xlabel('\epsilon_{zz} (10^{-3} /yr)'); title('Strain rate');
legend('inverted','truth','Location','SouthWest');

out_fn = fullfile(fig_dir,'synthetic_vertical_velocity.png');
print(h, out_fn, '-dpng', '-r120');
fprintf('Wrote %s\n', out_fn);
