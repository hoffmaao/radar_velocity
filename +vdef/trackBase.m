function base = trackBase(img, z, opts)
%TRACKBASE Follow the ice-base echo along a radargram.
%
%   base = vdef.trackBase(img, z, opts) returns the base depth for every
%   column of img (log power, depth rows z by along-track columns), as a
%   row vector in the units of z.
%
%   A windowed maximum fails where the base leaves the window - on the
%   EAGER lines it rises from ~285 m to ~210 m over the last 0.6 km, and a
%   250-350 m search snapped onto a flat system stripe at ~280 m instead.
%   So the base is TRACKED: start at the column opts.start (default the
%   middle), take the strongest return in opts.window there, then step
%   outward one column at a time, taking the strongest return within
%   +/- opts.step of the previous column's depth (and above opts.zmin). A
%   running median of opts.smooth columns removes single-column jumps, except
%   in the half-window at each end, where it would lag a sloping base.
%
%   opts.window [250 330], opts.step 6, opts.zmin 150, opts.smooth 9,
%   opts.start round(ncol/2); depth and step in the units of z.

if nargin < 3, opts = struct(); end
def = struct('window', [250 330], 'step', 6, 'zmin', 150, 'smooth', 9, 'start', []);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end
z = z(:); nc = size(img, 2);
sm = movmean(movmean(double(img), 5, 1, 'omitnan'), 3, 2, 'omitnan');
c0 = opts.start; if isempty(c0), c0 = round(nc/2); end
win = z >= opts.window(1) & z <= opts.window(2); zw = z(win);
[~, i0] = max(sm(win, c0)); base = nan(1, nc); base(c0) = zw(i0);
for dirn = [-1 1]
  c = c0 + dirn;
  while c >= 1 && c <= nc
    w = abs(z - base(c - dirn)) <= opts.step & z >= opts.zmin;
    zz = z(w); [~, j] = max(sm(w, c)); base(c) = zz(j);
    c = c + dirn;
  end
end
% running median for single-column jumps; at the line ends a median cannot be
% centred and shrinks onto one side, which lags a sloping base (the EAGER
% ramp runs to the end of the line), so the half-window at each end keeps
% the tracked values
raw = base; base = movmedian(raw, opts.smooth);
hw = min(floor(opts.smooth/2), floor(nc/2));
base(1:hw) = raw(1:hw); base(end-hw+1:end) = raw(end-hw+1:end);
end
