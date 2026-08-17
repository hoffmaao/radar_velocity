function s = apply_bulk_delay(s, tau, fc, fs)
%APPLY_BULK_DELAY Delay a slice by tau, envelope and carrier.
%   The same operation multipass.m:512-522 applies as z-motion
%   compensation, and the artefact vdef.coalignPair must undo. A delay tau
%   moves the envelope later by tau and multiplies by
%   exp(-1i*2*pi*fc*tau) in the matched-filter convention.
%
%   tau may be a scalar, or a 1 x Nx row giving a delay PER ALONG-TRACK
%   COLUMN. The per-column form is the realistic one: ref_z is a per-column
%   vector, so the misalignment multipass leaves behind varies along track.
%
%   Test support for test_vvel_task.m; a separate file rather than a
%   script-local function because Octave only registers script-local
%   functions once execution reaches their definition.
if isscalar(tau) && tau == 0
  return;
end
Nt = size(s, 1);
Nx = size(s, 2);
if isscalar(tau)
  tau = repmat(tau, 1, Nx);
else
  tau = tau(:).';
  assert(numel(tau) == Nx, 'tau has %d entries for %d columns', numel(tau), Nx);
end
dtb = 1/fs;
df = 1/(Nt*dtb);
f_bb = df * ifftshift(-floor(Nt/2):floor((Nt-1)/2)).';
S = fft(s, [], 1) .* exp(-1i*2*pi*(f_bb * tau));
s = ifft(S, [], 1);
s = bsxfun(@times, s, exp(-1i*2*pi*fc*tau));
end
