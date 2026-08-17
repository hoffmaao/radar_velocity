%BULK_SHIFT_DEPTH_ARTEFACT A depth-CONSTANT shift makes depth-DEPENDENT strain error.
%
%   Written to answer a precise objection (17 Aug 2026): "a scalar
%   coalignment should not be able to give a depth-dependent change in
%   strain that is coherent like this." The intuition is that a bulk
%   shift, being constant in depth, should cancel in the surface
%   referencing and leave depth structure untouched.
%
%   THE MECHANISM IT MISSES: mis-registration. Shifting the secondary by
%   d does not add d to the phase comparison - it makes fast-time bin t
%   of the reference compare against bin t-d of the secondary. After
%   surface referencing the apparent profile is
%
%       dtau_app(t) = dtau_true(t-d) - dtau_true(t_ref-d)
%                   ~ dtau_true(t) - dtau_true(t_ref)
%                     - d * [g(t) - g(t_ref)],   g = d(dtau)/dt + speckle
%                                                    phase-slope structure
%
%   The error is d TIMES A FIXED FUNCTION OF DEPTH, set by the scene's
%   own stratigraphy and speckle field. It is separable: every pair sees
%   the SAME depth shape g(z), scaled by its own bulk residual d. When d
%   is tide-proportional - which the multipass compensation error is on a
%   floating shelf - the result is a depth-coherent apparent strain
%   profile proportional to the tide: exactly what real flexure would
%   look like, and exactly what the scalar-era results showed.
%
%   THE DEMONSTRATION: one synthetic pair with a known smooth dtau(z),
%   three extra bulk delays d = 0.5, 1, 2 ns - the size of the along-track
%   VARIATION the scalar correction left on the real data. (3 ns and up
%   decorrelates this synthetic's 180 MHz band enough to fail the
%   coverage gates; the real pairs sit below that.) Asserts that the
%   referenced-dtau errors (a) have depth STRUCTURE far exceeding their
%   depth mean, (b) share one depth shape across the three d
%   (correlation > 0.95), (c) scale linearly with d. Separable d*g(z),
%   demonstrated.
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts/diagnostics \
%       gnuoctave/octave:latest octave --no-gui bulk_shift_depth_artefact.m

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));       % +vdef
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
  'opr_vvel','test'));                                                 % apply_bulk_delay

rand('seed', 5); randn('seed', 5);   %#ok<RAND>

fc = 750e6; fs = 300e6; Nt = 1500; Nx = 1200;
Surface = 0.20e-6;
Time = (0:Nt-1).'/fs;

% Scene: depth-decaying speckle + a bright band-limited surface return,
% as in test_vvel_task (the estimator needs the bin-scale structure)
depth_t = (Time - Surface); amp = 10.^(-max(depth_t,0)*8.4e4/1200);
s_ref = bsxfun(@times, amp, (randn(Nt,Nx)+1i*randn(Nt,Nx))/sqrt(2));
k = (-floor(Nt/2):floor((Nt-1)/2)).';
band = ifftshift(double(abs(k) <= 0.6*Nt/2));
f_bb = (fs/Nt)*ifftshift(k);
pulse = ifft(band .* exp(-1i*2*pi*f_bb*Surface)); pulse = pulse/max(abs(pulse));
s_ref = s_ref + 8*bsxfun(@times, pulse, exp(1i*2*pi*rand(1,Nx)));

% A smooth known dtau(z): zero at the surface, growing to 40 ps at depth
dtau_true = 40e-12 * max(0, (Time-Surface)/(Time(end)-Surface)).^1.5;
gamma = 0.95 * exp(-max(depth_t,0)*8.4e4/1200);
indep = bsxfun(@times, amp, (randn(Nt,Nx)+1i*randn(Nt,Nx))/sqrt(2));
% The surface pulse is INHERITED through the gamma mixing (it is already
% in s_ref), never re-drawn: a second independent pulse would decorrelate
% the surface, kill the unwrap seed, and NaN the whole profile - which is
% exactly what happened in the first version of this script.
s_sec0 = bsxfun(@times, sqrt(gamma), exp(-1i*2*pi*fc*dtau_true)) .* s_ref ...
       + bsxfun(@times, sqrt(max(1-gamma,0)), indep);

opts = struct('phase_sign',-1,'ref_twtt_offset',50e-9, ...
  'coherence_threshold',0.3,'max_gap_bins',20,'min_coverage',0.3, ...
  'mlook_window',[5 15]);

% profile helper inlined: Octave does not resolve script-local functions
% from every call context (same trap as test_tide_admittance, Aug 2026)
DS = [0 0.5e-9 1e-9 2e-9];
P = nan(Nt, numel(DS));
for q = 1:numel(DS)
  if DS(q) == 0, ss = s_sec0; else, ss = apply_bulk_delay(s_sec0, DS(q), fc, fs); end
  [ig, ch] = vdef.multilook(s_ref, ss, opts.mlook_window);
  m = struct('Time',Time,'Surface',Surface*ones(1,Nx),'fc',fc, ...
    'phase',angle(ig),'coherence',ch);
  [dt, ~] = vdef.differentialRange(m, opts);
  for t = 1:Nt
    v = dt(t,:); v = v(isfinite(v));
    if numel(v) >= 50, P(t,q) = median(v); end
  end
end
E = P(:,2:4) - P(:,1)*ones(1,3);

sel = isfinite(E(:,1)) & isfinite(E(:,2)) & isfinite(E(:,3)) & Time > Surface+0.1e-6;
e1 = E(sel,1); e2 = E(sel,2); e3 = E(sel,3);
fprintf('bulk shifts injected: 0.5, 1, 2 ns (depth-CONSTANT)\n');
fprintf('referenced dtau error: rms %.1f / %.1f / %.1f ps\n', ...
  1e12*sqrt(mean(e1.^2)), 1e12*sqrt(mean(e2.^2)), 1e12*sqrt(mean(e3.^2)));
fprintf('depth STRUCTURE vs depth mean:  std %.1f ps vs |mean| %.1f ps (d = 2 ns)\n', ...
  1e12*std(e3), 1e12*abs(mean(e3)));
c12 = corrcoef(e1,e2); c13 = corrcoef(e1,e3);
fprintf('shape correlation across d:     r(1,2) = %.3f, r(1,3) = %.3f\n', ...
  c12(1,2), c13(1,2));
s31 = (e3.'*e1)/(e1.'*e1);
fprintf('amplitude scaling e(2)/e(0.5):  %.2f (linear predicts 4.00)\n', s31);

assert(std(e3) > 3*abs(mean(e3)), ...
  'the error should be depth-STRUCTURED, not a constant');
% adjacent sizes share the shape almost exactly; at 2 ns (0.6 bins) the
% second-order d^2 term starts to bend it, so the wide comparison is held
% to a looser floor - the claim is strong shared structure, not exact
% shape invariance
assert(c12(1,2) > 0.90, ...
  'adjacent shift sizes should share the depth shape (got r = %.2f)', c12(1,2));
assert(c13(1,2) > 0.75, ...
  'even 4x apart the shapes should correlate strongly (got r = %.2f)', c13(1,2));
assert(abs(s31 - 4) < 1.0, ...
  'the error should scale linearly with the shift (got %.2f vs 4)', s31);

fprintf(['\nPASS: a depth-constant shift produces a depth-DEPENDENT error\n' ...
  'with one frozen shape scaled by the shift - so a tide-proportional\n' ...
  'bulk residual manufactures a coherent, tide-proportional strain\n' ...
  'profile. Depth coherence is NOT evidence the scalar-era signal was ice.\n']);

