function [igram, coh] = multilook(s_ref, s_sec, win)
%MULTILOOK Boxcar-multilooked repeat-pass interferogram and coherence.
%   [igram, coh] = MULTILOOK(s_ref, s_sec, win) forms
%
%     igram = < s_sec .* conj(s_ref) >
%     coh   = |< s_sec conj(s_ref) >| / sqrt(<|s_ref|^2> <|s_sec|^2>)
%
%   where <.> is a boxcar average over win = [Nt_win Nx_win] samples
%   (fast time x along track). Both inputs are Nt x Nx complex SLC images
%   already coregistered onto a common time and along-track axis - i.e.
%   the pass slices of the CSARP_multipass comp_mode 3 product.
%
%   ORDERING RULE (carried over from the delta-k work in the fabric
%   project): the cross product is formed PER PIXEL and only then averaged.
%   Averaging the individual images first would decorrelate exactly where
%   the fringe rate is high.
%
%   Non-finite samples are excluded from each window rather than poisoning
%   it, so edges and data gaps degrade gracefully.
%
%   See also vdef.differentialRange.

if nargin < 3 || isempty(win)
  win = [5 15];
end

s_ref = double(s_ref);
s_sec = double(s_sec);

good = isfinite(s_ref) & isfinite(s_sec);
s_ref(~good) = 0;
s_sec(~good) = 0;

cross = s_sec .* conj(s_ref);

num   = boxsum(cross,        win);
p_ref = boxsum(abs(s_ref).^2, win);
p_sec = boxsum(abs(s_sec).^2, win);
cnt   = boxsum(double(good), win);

cnt(cnt < 1) = NaN;
igram = num ./ cnt;

den = sqrt(p_ref .* p_sec);
coh = abs(num) ./ den;
coh(~isfinite(coh)) = 0;
coh = min(coh, 1);          % guard against round-off above unity

end

function y = boxsum(x, win)
% Separable boxcar SUM with 'same' support (zeros outside), NaN-free input.
kt = ones(max(1,round(win(1))), 1);
kx = ones(1, max(1,round(win(2))));
y  = conv2(kt, kx, x, 'same');
end
