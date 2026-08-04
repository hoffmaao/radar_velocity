function [dtau, info] = differentialRange(map, opts)
%DIFFERENTIALRANGE Differential two-way traveltime from a repeat-pass pair.
%   [dtau, info] = DIFFERENTIALRANGE(map, opts) converts the interferogram
%   phase of a repeat-pass pair into the change in two-way traveltime
%   dtau(twtt, x) = tau_sec - tau_ref, referenced to zero just below the
%   surface return.
%
%   map fields:
%     .Time        Nt x 1 fast-time axis of the coregistered pair [s]
%     .Surface     1 x Nx surface twtt [s]
%     .fc          centre frequency [Hz]
%     .phase       Nt x Nx interferogram phase [rad] (wrapped, or already
%                  unwrapped if .phase_is_unwrapped is true)
%     .coherence   Nt x Nx magnitude coherence in [0, 1]
%     .phase_is_unwrapped  logical (default false)
%
%   opts fields (see vdef.m for the full documented defaults):
%     .phase_sign            +1 or -1; dtau = phase_sign*Phi/(2*pi*fc)
%     .ref_twtt_offset       twtt below the surface where dtau := 0 [s]
%     .coherence_threshold   samples below this do not constrain the unwrap
%     .max_gap_bins          longest incoherent run that may be bridged
%
%   SIGN CONVENTION. The interferogram is formed as sec .* conj(ref) (see
%   vdef.multilook), matching both polarimetric.m in OPR and the
%   fabric_anisotropy chain. For the matched-filter convention
%   s ~ exp(-1i*2*pi*fc*tau) this gives Phi = -2*pi*fc*dtau, hence
%   phase_sign = -1, which is the value settled on in the fabric project by
%   regressing phase against coregistration offsets. There is no
%   coregistration field in the multipass product to re-derive it from, so
%   it is a fixed parameter here; vdef.m defaults it to -1 and the
%   validation path is the cross-check described in opr_vvel/README.md.
%
%   UNWRAPPING. dtau is smooth in depth, so the phase is unwrapped along
%   fast time outward from the surface reference bin rather than with a 2-D
%   unwrapper. Wrapped fast-time differences are accumulated; steps whose
%   endpoints are below the coherence threshold contribute zero (i.e. the
%   gap is assumed not to cross a fringe), and once a contiguous incoherent
%   run exceeds max_gap_bins everything deeper in that trace is marked
%   invalid rather than silently extrapolated. This is deliberately
%   conservative: it fails toward "no data" instead of accumulating
%   confident nonsense the way a region-growing unwrapper does.
%
%   Returns:
%     dtau       Nt x Nx differential traveltime [s], NaN where invalid
%     info.ref_bin        1 x Nx index of the surface reference bin
%     info.valid          Nt x Nx logical validity mask
%     info.max_valid_bin  1 x Nx deepest trusted bin per trace
%     info.phase_sign     the sign actually applied
%
%   See also vdef.multilook, vdef.blockAverage, vdef.verticalDisplacement.

Time = map.Time(:);
Nt   = numel(Time);
Nx   = size(map.phase, 2);
dt   = Time(2) - Time(1);

phase_is_unwrapped = isfield(map,'phase_is_unwrapped') && map.phase_is_unwrapped;

coh = map.coherence;
if isempty(coh)
  coh = ones(Nt, Nx);
end
coherent = coh >= opts.coherence_threshold;

%% Surface reference bin per trace
Surface = map.Surface(:).';
ref_twtt = Surface + opts.ref_twtt_offset;
ref_bin  = round(interp1(Time, 1:Nt, ref_twtt, 'linear', NaN));
ref_bin  = min(max(ref_bin, 1), Nt);

%% Unwrap along fast time, outward from the reference bin
valid = true(Nt, Nx);
max_valid_bin = repmat(Nt, 1, Nx);

if phase_is_unwrapped
  Phi = map.phase;
  valid = coherent;
else
  Phi = nan(Nt, Nx);
  % Wrapped first difference along fast time
  dphi = angle(exp(1i*(map.phase(2:end,:) - map.phase(1:end-1,:))));
  % A step is only trusted when BOTH of its endpoints are coherent
  step_ok = coherent(2:end,:) & coherent(1:end-1,:);
  dphi(~step_ok) = 0;

  for x = 1:Nx
    r = ref_bin(x);
    if ~isfinite(r)
      valid(:,x) = false;
      max_valid_bin(x) = 0;
      continue;
    end
    Phi(r,x) = map.phase(r,x);

    % Downward (increasing twtt)
    if r < Nt
      Phi(r+1:Nt,x) = map.phase(r,x) + cumsum(dphi(r:Nt-1,x));
    end
    % Upward (decreasing twtt); only used for QC, not for the inversion
    if r > 1
      Phi(1:r-1,x) = map.phase(r,x) - flipud(cumsum(flipud(dphi(1:r-1,x))));
    end

    % Truncate below the first over-long incoherent run
    stop = firstLongGap(step_ok(r:Nt-1,x), opts.max_gap_bins);
    if ~isempty(stop)
      last = r + stop - 1;
      valid(last+1:Nt,x) = false;
      max_valid_bin(x)   = last;
    end
  end
  valid = valid & coherent;
end

%% Reference to zero at the surface bin and convert to traveltime
idx = ref_bin + (0:Nx-1)*Nt;
ok  = isfinite(idx);
Phi_ref = nan(1, Nx);
Phi_ref(ok) = Phi(idx(ok));
Phi = bsxfun(@minus, Phi, Phi_ref);

dtau = opts.phase_sign * Phi / (2*pi*map.fc);
dtau(~valid) = NaN;

% Nothing above the reference bin is used: there is no differential column
% between the surface and a reflector shallower than the reference.
for x = 1:Nx
  if isfinite(ref_bin(x)) && ref_bin(x) > 1
    dtau(1:ref_bin(x)-1, x) = NaN;
  end
end

info = [];
info.ref_bin       = ref_bin;
info.valid         = isfinite(dtau);
info.max_valid_bin = max_valid_bin;
info.phase_sign    = opts.phase_sign;
info.dt            = dt;

end

function stop = firstLongGap(step_ok, max_gap_bins)
% Index (into step_ok) of the start of the first run of >max_gap_bins
% consecutive untrusted steps; empty if there is none.
stop = [];
if isempty(step_ok) || ~isfinite(max_gap_bins) || max_gap_bins < 0
  return;
end
bad = ~step_ok(:).';
if ~any(bad)
  return;
end
d  = diff([0 bad 0]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
runlen = stops - starts + 1;
first  = find(runlen > max_gap_bins, 1, 'first');
if ~isempty(first)
  stop = starts(first);
end
end
