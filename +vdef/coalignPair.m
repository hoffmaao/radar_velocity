function [s_sec, info] = coalignPair(s_ref, s_sec, map, opts)
%COALIGNPAIR Measure and remove the bulk fast-time shift between a pair.
%   [s_sec, info] = COALIGNPAIR(s_ref, s_sec, map, opts) measures the bulk
%   misalignment of the secondary slice relative to the reference slice in
%   a window around the surface return, to well under a bin, and applies
%   the inverse time shift (envelope and carrier phase) to the secondary
%   so the pair is internally aligned before the interferogram is formed.
%
%   WHY THIS EXISTS. multipass comp_mode 3 motion-compensates each pass's
%   FCS z-motion by a time shift of ref_z/(c/2) (multipass.m:512-522,
%   envelope and phase), and pass.surface is not updated for that shift
%   (layers.twtt_ref explicitly is, multipass.m:667-669). On grounded ice
%   ref_z is real platform motion and the compensation is correct. On a
%   FLOATING shelf the platform and the surface ride the tide together,
%   the antenna-to-surface range does not change, and ref_z is essentially
%   the tide - so the compensation displaces the returns by a
%   tide-proportional amount that never was a range change. Measured with
%   this estimator on the EAGER 2022 products, the residual misalignment
%   reaches 7 ns (about 2 range bins) and runs at 1.2-1.3 times
%   -(ref_z_sec - ref_z_ref)/(c/2), i.e. slightly MORE than the full
%   erroneous compensation survives to the product. It contaminated the
%   inferred strain at roughly 57 mm of apparent column displacement per
%   metre of tide before this fix.
%
%   The correction is EMPIRICAL, from the data itself, rather than the
%   deterministic inverse: the surviving fraction is not exactly one and
%   is not knowable in advance, and an empirical measure also absorbs any
%   per-pair coregistration_time_shift applied upstream.
%
%   THE ESTIMATOR is the group delay from the phase slope of the
%   cross-spectrum of the two windows: for b(t) = a(t - tau),
%   P(f) = B(f) conj(A(f)) = |A(f)|^2 exp(-1i*2*pi*f*tau), so consecutive
%   frequency bins step in phase by -2*pi*df*tau. Summing
%   P(f+df) conj(P(f)) over the band and over traces gives tau from a
%   single angle, unambiguous over +/- W*dt/2 (window length W). This is
%   deliberately NOT a parabolic refinement of the envelope correlation
%   peak: on speckle the correlation peak is about one bin wide, and a
%   parabola through three integer-lag samples of a one-bin peak
%   underestimates sub-bin shifts several-fold. The carrier phase
%   exp(-1i*2*pi*fc*tau) is constant across the band and cancels in the
%   slope, as does any spectral shape common to the two slices.
%
%   map fields: .Time (Nt x 1) [s], .Surface (1 x Nx) [s], .fc [Hz]
%   opts fields (defaults in opr_vvel/vvel_defaults.m):
%     .coalign_max_lag      largest credible shift [bins]; a larger
%                           estimate is rejected as a failed measurement
%     .coalign_half_win     half-width of the surface window [bins]
%     .coalign_min_quality  minimum info.quality for the estimate to be
%                           applied; below it the measurement is rejected
%                           as noise and the pair is left unaligned
%
%   Returns:
%     s_sec            the secondary, advanced by the measured delay
%     info.dtau_bulk   measured extra delay of sec relative to ref [s]
%     info.quality     phase-slope consistency of the cross-spectrum, in
%                      [0..1]; low values mean the estimate is noise
%     info.applied     false when no reliable estimate was possible; the
%                      secondary is then returned unchanged and dtau_bulk
%                      is NaN. A quality-floor rejection still reports the
%                      measured quality so the product records why.
%
%   SIGN. Positive dtau_bulk means the secondary's returns arrive LATER
%   than the reference's. In the matched-filter convention a delay tau
%   multiplies the signal by exp(-1i*2*pi*fc*tau), so the correction
%   multiplies by exp(+1i*2*pi*fc*dtau_bulk) and advances the envelope by
%   dtau_bulk.
%
%   See also vdef.multilook, vdef.differentialRange.

info = struct('dtau_bulk', NaN, 'quality', NaN, 'applied', false);

Time = map.Time(:);
Nt = size(s_ref, 1);
dt = Time(2) - Time(1);
maxlag  = max(1, round(opts.coalign_max_lag));
halfwin = max(4, round(opts.coalign_half_win));
minq    = opts.coalign_min_quality;

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
win = max(1, b0-halfwin) : min(Nt, b0+halfwin);
W = numel(win);
if W < 16
  warning('coalignPair: surface window of %d bins is too short; pair left unaligned.', W);
  return;
end

% Complex fields in the surface window, tapered so the window edges do not
% leak into the cross-spectrum phase; non-finite samples contribute zero
a = double(s_ref(win,:));
b = double(s_sec(win,:));
a(~isfinite(a)) = 0;
b(~isfinite(b)) = 0;
taper = 0.5 - 0.5*cos(2*pi*(0:W-1).'/(W-1));   % hann, no toolbox dependency
Fa = fft(bsxfun(@times, a, taper));
Fb = fft(bsxfun(@times, b, taper));
P = sum(Fb .* conj(Fa), 2);                    % cross-spectrum over traces

if ~any(isfinite(P)) || ~any(abs(P) > 0)
  warning('coalignPair: empty surface window; pair left unaligned.');
  return;
end

% Group delay from the phase step between consecutive frequency bins.
% Pairs near the Nyquist wrap (the middle of fft order) are EXCLUDED: the
% delay was applied over the whole record, not cyclically over the window,
% so the phase ramp is discontinuous across the band edge, and window
% leakage concentrates enough magnitude there to drag the vector average
% several-fold (verified: without the guard an injected 1.30 ns measured
% as 0.27 ns; with it, 1.299-1.300 ns across the tested range).
pp = P(2:end) .* conj(P(1:end-1));
kpair = (1:W-1).';
keep = abs(kpair - W/2) > 0.15*W;
q = sum(pp(keep));
dfw = 1/(W*dt);
dtau = -angle(q) / (2*pi*dfw);
quality = abs(q) / max(sum(abs(pp(keep))), eps);

if ~isfinite(dtau)
  warning('coalignPair: group delay estimate is not finite; pair left unaligned.');
  return;
end
if quality < minq
  warning('coalignPair: cross-spectrum quality %.2f is below the %.2f floor; pair left unaligned.', ...
    quality, minq);
  info.quality = quality;
  return;
end
if abs(dtau) > maxlag*dt
  warning('coalignPair: measured %.2f ns exceeds the %.2f ns credibility bound; pair left unaligned.', ...
    dtau*1e9, maxlag*dt*1e9);
  return;
end

% Advance the secondary by the measured delay: envelope through the
% baseband spectrum, carrier phase at fc
df = 1/(Nt*dt);
f_bb = df * ifftshift(-floor(Nt/2):floor((Nt-1)/2)).';
s_sec = ifft(bsxfun(@times, fft(double(s_sec)), exp(1i*2*pi*f_bb*dtau))) ...
  * exp(1i*2*pi*map.fc*dtau);

info.dtau_bulk = dtau;
info.quality   = quality;
info.applied   = true;

end
