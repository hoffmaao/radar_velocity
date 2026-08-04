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

tol = 0.20;
assert(abs(mean(S1_err,'omitnan')) < tol*abs(S1_true), ...
  'S1 not recovered within %.0f%%', 100*tol);
assert(abs(mean(S2_err,'omitnan')) < tol*abs(S2_true), ...
  'S2 not recovered within %.0f%%', 100*tol);
fprintf('\nPASS: strain rate recovered within %.0f%% of truth.\n', 100*tol);

%% 5. Figure
% =====================================================================
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

h = figure('Visible','off');

subplot(1,3,1);
plot(1e12*blk.dtau, V.depth(:,1), '-'); hold on;
plot(1e12*dtau_true, depth_col, 'k--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); ylim([0 550]); grid on;
xlabel('\Delta\tau (ps)'); ylabel('Depth (m)');
title('Differential traveltime');
legend('blocks','truth','Location','SouthWest');

subplot(1,3,2);
plot(V.v, V.depth(:,1), '-'); hold on;
plot(S.v_fit, S.depth_grid, 'k-', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); ylim([0 550]); grid on;
xlabel('Relative w (m/yr)'); title('Vertical velocity');

subplot(1,3,3);
plot(S.eps_zz, S.depth_grid, '-', 'LineWidth', 1.5); hold on;
dd = linspace(0,H_norm,100);
plot(eps_zz_true(dd), dd, 'k--', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); ylim([0 550]); grid on;
xlabel('\epsilon_{zz} (1/yr)'); title('Strain rate');
legend('inverted','truth','Location','SouthWest');

out_fn = fullfile(fig_dir,'synthetic_vertical_velocity.png');
print(h, out_fn, '-dpng', '-r120');
fprintf('Wrote %s\n', out_fn);
