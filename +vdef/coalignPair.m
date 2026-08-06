function [s_sec, info] = coalignPair(s_ref, s_sec, map, opts)
%COALIGNPAIR Measure and remove the fast-time misalignment between a pair.
%   [s_sec, info] = COALIGNPAIR(s_ref, s_sec, map, opts) measures how far
%   the secondary slice is displaced from the reference in fast time, as a
%   function of ALONG-TRACK POSITION, and applies the inverse shift per
%   column (envelope and carrier phase) so the pair is internally aligned
%   before the interferogram is formed.
%
%   WHY THIS EXISTS. multipass comp_mode 3 motion-compensates each pass's
%   FCS z-motion by a time shift of ref_z/(c/2) (multipass.m:512-522,
%   envelope and phase), and pass.surface is not updated for that shift
%   (layers.twtt_ref explicitly is, multipass.m:667-669). On grounded ice
%   ref_z is real platform motion and the compensation is correct. On a
%   FLOATING shelf the platform and the surface ride the tide together,
%   the antenna-to-surface range does not change, and ref_z is essentially
%   the tide - so the compensation displaces the returns by a
%   tide-proportional amount that never was a range change.
%
%   WHY IT MUST BE PER COLUMN. ref_z is a per-column vector, not a number.
%   The two passes walk the line at different times and sometimes in
%   opposite directions, so the tide changes by different amounts at each
%   end and the misalignment varies ALONG TRACK. On the EAGER pairs its
%   5th-to-95th-percentile spread within a single pair reaches 1.1 m -
%   7.4 ns, over two range bins - which is LARGER than the mean offset.
%
%   A scalar correction is not merely incomplete here, it is actively
%   misleading. Both this function's earlier scalar form and multipass's
%   own coregistration_time_shift remove a line MEAN, so what survives is
%   proportional to (a(x) - 1), where a(x) is the local surface tidal
%   admittance normalised to line-mean 1. That residual changes sign
%   exactly where a(x) = 1 - a position fixed by the arbitrary
%   normalisation of the survey, not by the ice - and it duly manufactured
%   an apparent flexure hinge within 0.5 km of that crossing in all five
%   EAGER products. Removing the along-track structure is what makes the
%   inferred strain mean anything.
%
%   WHETHER ANY OF THIS REACHES A GIVEN PRODUCT DEPENDS ON ITS
%   CALIBRATION. param.multipass.coregistration_time_shift applies a
%   per-pass fast-time shift from the comp_mode 2 stage and removes the
%   MEAN of this misalignment as a side effect. EAGER_2022 and GL4 carry
%   all zeros and show the full offset (0.97-1.07 times the prediction);
%   GL1, GL2 and GL3 carry nonzero shifts and are aligned in the mean.
%   None of them are aligned in the along-track VARIATION, because
%   coregistration_time_shift is also only a scalar.
%
%   The correction is EMPIRICAL, measured from the data rather than
%   computed from ref_z, so it absorbs whatever coregistration_time_shift
%   already did and does not depend on the compensation surviving at any
%   particular fraction.
%
%   THE ESTIMATOR, per along-track window, is the normalised
%   cross-correlation of the two slices' TRACE-AVERAGED POWER profiles
%   over a window around the surface, FFT-upsampled by
%   opts.coalign_upsample so the peak is located to a fraction of a bin
%   without a parabolic fit.
%
%   Why the trace average is well conditioned: a single trace's envelope
%   is speckle and its correlation peak is about one bin wide, which is
%   what made a parabolic peak fit underestimate sub-bin shifts several
%   fold. Averaging |s|^2 over the columns in the window first leaves the
%   smooth pulse shape of the surface return. The surface twtt varies
%   along track, which smears that average - but both slices are resampled
%   onto the SAME along-track axis by multipass, so the smearing is common
%   to the two profiles and cancels in their cross-correlation.
%
%   Why the envelope rather than the carrier: at fc = 750 MHz one carrier
%   cycle is 1.33 ns, so the 5 ns shifts seen here span about four cycles
%   and the carrier phase alone cannot resolve them. The misalignment is a
%   genuine full time shift, so the envelope determines it, and the
%   correction is then applied to envelope and carrier together.
%
%   THE LAG BOUND IS LOAD-BEARING, not a formality. The trace-averaged
%   surface profile carries a range sidelobe of the transmit pulse at
%   about +/-5.4 bins, at 0.37-0.55 of the main peak in every EAGER pair.
%   On one pair (GL1 / 20221211_07) that sidelobe outranked the true peak
%   and produced a -19 ns estimate against a true 0.8 ns. A peak-dominance
%   test does NOT separate these cases: that pair's peak-to-sidelobe ratio
%   was 1.8 while a correctly measured pair (GL1 / 20221211_01) sat at
%   1.03. What does separate them is physics - the misalignment cannot
%   exceed the tidal range over c/2 (about 3 bins here) - so
%   opts.coalign_max_lag defaults to 3 bins, comfortably above every true
%   shift measured and safely inside the sidelobe. A peak found AT the
%   search boundary is rejected rather than clamped.
%
%   map fields: .Time (Nt x 1) [s], .Surface (1 x Nx) [s], .fc [Hz]
%   opts fields (defaults in opr_vvel/vvel_defaults.m):
%     .coalign_max_lag      largest credible shift [bins]; a peak at or
%                           beyond this is rejected as a failed measurement
%     .coalign_half_win     half-width of the surface window [bins]
%     .coalign_upsample     cross-correlation upsampling factor
%     .coalign_min_quality  minimum window correlation for a window's
%                           estimate to be used
%     .coalign_win_cols     along-track window length [columns]. Empty, 0
%                           or Inf gives ONE window over the whole line,
%                           i.e. the old scalar behaviour, which is kept
%                           only for reproducing the artefact deliberately.
%
%   Returns:
%     s_sec              the secondary, advanced by the measured delay
%     info.dtau_bulk     mean of the applied profile [s] - the scalar that
%                        the earlier form of this function reported, kept
%                        so existing products and analyses stay readable
%     info.dtau_profile  1 x Nx applied delay [s]
%     info.dtau_win      per-window measured delay [s]
%     info.x_win         per-window centre column
%     info.quality_win   per-window normalised correlation
%     info.quality       median of the accepted windows' correlation
%     info.peak_ratio    worst (largest) window sidelobe-to-peak ratio, for
%                        QC only - see above for why it must not gate
%     info.n_win, info.n_win_ok   windows attempted and accepted
%     info.applied       false when no reliable estimate was possible; the
%                        secondary is then returned unchanged. A
%                        quality-floor rejection still reports the measured
%                        quality so the product records why.
%
%   SIGN. Positive delay means the secondary's returns arrive LATER than
%   the reference's. In the matched-filter convention a delay tau
%   multiplies the signal by exp(-1i*2*pi*fc*tau), so the correction
%   multiplies by exp(+1i*2*pi*fc*tau) and advances the envelope by tau.
%
%   See also vdef.multilook, vdef.differentialRange.

Nt = size(s_ref, 1);
Nx = size(s_ref, 2);

info = struct('dtau_bulk', NaN, 'dtau_profile', [], 'dtau_win', [], ...
  'x_win', [], 'quality_win', [], 'quality', NaN, 'peak_ratio', NaN, ...
  'n_win', 0, 'n_win_ok', 0, 'applied', false);

Time = map.Time(:);
dt = Time(2) - Time(1);
maxlag  = max(1, opts.coalign_max_lag);
halfwin = max(4, round(opts.coalign_half_win));
minq    = opts.coalign_min_quality;
if isfield(opts,'coalign_upsample') && ~isempty(opts.coalign_upsample)
  UP = max(1, round(opts.coalign_upsample));
else
  UP = 32;
end
if isfield(opts,'coalign_win_cols') && ~isempty(opts.coalign_win_cols) ...
    && isfinite(opts.coalign_win_cols) && opts.coalign_win_cols > 0
  wcols = min(Nx, round(opts.coalign_win_cols));
else
  wcols = Nx;                     % one window: the legacy scalar behaviour
end

% Fast-time window around the surface
sfc = mean(map.Surface, 'omitnan');
if ~isfinite(sfc)
  warning('coalignPair: no finite surface twtt; pair left unaligned.');
  return;
end
b0 = round(interp1(Time, 1:Nt, sfc, 'linear', NaN));
if ~isfinite(b0)
  warning('coalignPair: surface twtt outside the fast-time axis; pair left unaligned.');
  return;
end
tw = max(1, b0-halfwin) : min(Nt, b0+halfwin);
W = numel(tw);
if W < 16
  warning('coalignPair: surface window of %d bins is too short; pair left unaligned.', W);
  return;
end
taper = 0.5 - 0.5*cos(2*pi*(0:W-1).'/(W-1));   % hann, no toolbox dependency

% Along-track windows, 50% overlap. Each is averaged over wcols columns, so
% the estimates are already smooth; they are interpolated, not filtered.
step   = max(1, round(wcols/2));
starts = 1:step:max(1, Nx-wcols+1);
if isempty(starts), starts = 1; end
if starts(end) + wcols - 1 < Nx && wcols < Nx
  starts(end+1) = Nx - wcols + 1;      % cover the tail
end
nw = numel(starts);

dtau_win = nan(1, nw); q_win = nan(1, nw); pr_win = nan(1, nw); x_win = nan(1, nw);
for k = 1:nw
  cols = starts(k) : min(Nx, starts(k)+wcols-1);
  x_win(k) = mean(cols);
  p_ref = trace_mean_power(s_ref(tw, cols), taper);
  p_sec = trace_mean_power(s_sec(tw, cols), taper);
  [db, qq, pr] = window_shift(p_ref, p_sec, UP, maxlag);
  if ~isfinite(db), continue; end
  pr_win(k) = pr;
  q_win(k)  = qq;
  if qq < minq
    continue;                     % measured, but not trusted
  end
  if abs(db) >= maxlag - 1/UP
    continue;                     % peak on the credibility bound
  end
  dtau_win(k) = db * dt;
end

info.n_win       = nw;
info.x_win       = x_win;
info.dtau_win    = dtau_win;
info.quality_win = q_win;
info.peak_ratio  = max_or_nan(pr_win);

ok = isfinite(dtau_win);
info.n_win_ok = nnz(ok);
if ~any(ok)
  qmax = max_or_nan(q_win);
  info.quality = qmax;
  warning('coalignPair: no along-track window produced a usable estimate (best correlation %.2f against a %.2f floor); pair left unaligned.', ...
    qmax, minq);
  return;
end
info.quality = median(q_win(ok));

% A single usable window cannot describe along-track structure; treat its
% estimate as a constant rather than pretending to a profile.
if nnz(ok) == 1
  dtau_x = repmat(dtau_win(ok), 1, Nx);
  if nw > 1
    warning('coalignPair: only 1 of %d along-track windows was usable, so a CONSTANT shift is being applied and the along-track residual is NOT corrected.', nw);
  end
else
  % Linear between window centres, held flat beyond the outermost ones -
  % extrapolating a slope past the data would invent shifts at the ends
  xv = x_win(ok); dv = dtau_win(ok);
  dtau_x = interp1(xv, dv, 1:Nx, 'linear');
  dtau_x(1:floor(xv(1)))  = dv(1);
  dtau_x(ceil(xv(end)):Nx) = dv(end);
  dtau_x(~isfinite(dtau_x)) = 0;
end

% Advance the secondary by the measured delay: envelope through the
% baseband spectrum, carrier phase at fc, both per column
df   = 1/(Nt*dt);
f_bb = df * ifftshift(-floor(Nt/2):floor((Nt-1)/2)).';
S = fft(double(s_sec), [], 1);
S = S .* exp(1i*2*pi*(f_bb * dtau_x));           % Nt x 1 times 1 x Nx
s_sec = ifft(S, [], 1);
clear S;
s_sec = bsxfun(@times, s_sec, exp(1i*2*pi*map.fc*dtau_x));

info.dtau_profile = dtau_x;
info.dtau_bulk    = mean(dtau_x, 'omitnan');
info.applied      = true;

end

%% ========================================================================
function [dtau_bins, peak, pratio] = window_shift(p_ref, p_sec, UP, maxlag)
% Sub-bin shift of p_sec relative to p_ref from the normalised, FFT-
% upsampled cross-correlation, searched only within +/- maxlag bins.
dtau_bins = NaN; peak = NaN; pratio = NaN;

den = sqrt(sum(p_ref.^2) * sum(p_sec.^2));
if ~(den > 0) || ~all(isfinite([p_ref; p_sec]))
  return;
end

W    = numel(p_ref);
Nfft = 2^nextpow2(2*W);
X  = fft(p_sec, Nfft) .* conj(fft(p_ref, Nfft));
Xu = zeros(Nfft*UP, 1);
h  = Nfft/2;
Xu(1:h)         = X(1:h);
Xu(end-h+1:end) = X(h+1:end);
xc = real(ifft(Xu)) * UP / den;

lags = ifftshift((-Nfft*UP/2 : Nfft*UP/2-1).') / UP;
[lags, ord] = sort(lags);
xc = xc(ord);

sel = abs(lags) <= maxlag;
if ~any(sel), return; end
lag_in = lags(sel);
[peak, im] = max(xc(sel));
dtau_bins = lag_in(im);
if ~isfinite(peak) || ~isfinite(dtau_bins)
  dtau_bins = NaN; peak = NaN; return;
end

% Sidelobe measured WIDER than the search bound: the sidelobe worth
% knowing about sits at about +/-5.4 bins, outside the 3-bin bound that
% deliberately excludes it, so a ratio computed only inside the bound
% would almost always come back NaN and tell nobody anything.
qc = abs(lags) <= max(2*maxlag, 8);
pratio = sidelobe_ratio(xc(qc), lags(qc), dtau_bins, peak);
end

%% ========================================================================
function p = trace_mean_power(s, taper)
% Mean power over the along-track columns, mean-removed and tapered.
% Non-finite samples are ignored rather than zeroed, so a column with a
% dropout does not pull the profile toward zero at that bin.
p = mean(abs(double(s)).^2, 2, 'omitnan');
p(~isfinite(p)) = 0;
p = (p - mean(p)) .* taper;
end

%% ========================================================================
function r = sidelobe_ratio(xc, lags, peak_lag, peak)
% Height of the tallest local maximum that is not the main peak, relative
% to the peak. Local maxima within half a bin of the peak lag are part of
% the main lobe, not competitors, so they are excluded by position.
r = NaN;
if ~(peak > 0), return; end
isloc = false(size(xc));
isloc(2:end-1) = xc(2:end-1) > xc(1:end-2) & xc(2:end-1) > xc(3:end);
isloc = isloc & abs(lags - peak_lag) > 0.5;
if ~any(isloc), return; end
r = max(xc(isloc)) / peak;
end

%% ========================================================================
function v = max_or_nan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end
