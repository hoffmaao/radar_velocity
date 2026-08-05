function s = apply_bulk_delay(s, tau, fc, fs)
%APPLY_BULK_DELAY Delay a whole slice by tau, envelope and carrier.
%   The same operation multipass.m:512-522 applies as z-motion
%   compensation, and the artefact vdef.coalignPair must undo. A delay tau
%   moves the envelope later by tau and multiplies by
%   exp(-1i*2*pi*fc*tau) in the matched-filter convention.
%
%   Test support for test_vvel_task.m; a separate file rather than a
%   script-local function because Octave only registers script-local
%   functions once execution reaches their definition.
if tau == 0
  return;
end
Nt = size(s, 1);
dtb = 1/fs;
df = 1/(Nt*dtb);
f_bb = df * ifftshift(-floor(Nt/2):floor((Nt-1)/2)).';
s = ifft(bsxfun(@times, fft(s), exp(-1i*2*pi*f_bb*tau))) * exp(-1i*2*pi*fc*tau);
end
