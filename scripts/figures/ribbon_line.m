function nsig = ribbon_line(ax, x, y, v, sg, dmap, CLIM, NSIG, mid, lw)
%RIBBON_LINE Colour-coded polyline whose SATURATION carries significance.
%
%   nsig = RIBBON_LINE(ax, x, y, v, sg, dmap, CLIM, NSIG, mid, lw) draws the
%   polyline (x, y) one segment at a time, colouring each segment by the
%   value v through the diverging map dmap over [-CLIM, +CLIM], and returns
%   how many segments reached NSIG sigma.
%
%   (x, y) is where the line goes; v is what it means. On the map those are
%   different things - the line follows a leg's track and v is the tidal
%   response along it. On a section plot y and v are the same array, and the
%   colour is then redundant with height BY DESIGN: it is what keys the
%   section to the map, and it carries the significance, which height alone
%   cannot show.
%
%   WHY SATURATION AND NOT COLOUR ALONE. Per-leg sigma reaches 3.3 mm/m
%   against a +/-4 colour scale, so a one-sigma excursion saturated the ramp
%   and a noisy leg was indistinguishable from a real signal. Each segment is
%   blended toward the neutral midpoint `mid` by w = min(1, |v|/(NSIG*sg)),
%   so a sample the data cannot separate from zero is drawn as very nearly
%   zero. Full colour means at least NSIG sigma.
%
%   THE COST, stated because it is real: a genuine zero and an unresolved
%   sample now look alike. That is the honest ambiguity - neither is evidence
%   of a response - but "pale" does not distinguish "no signal" from "cannot
%   tell", and a caption should not imply it does.
%
%   Segments touching a gap are skipped, so the line breaks over a stretch
%   with no data rather than bridging it. Colour is looked up through the
%   SAME index arithmetic the colour bar uses, so a ribbon and the bar cannot
%   drift apart.
%
%   See also tidal_response_map, line_sections.

ndiv = size(dmap, 1);
nsig = 0;
for k = 1:numel(x)-1
  if ~all(isfinite([x(k) x(k+1) y(k) y(k+1) v(k) v(k+1)])), continue; end
  vm = (v(k) + v(k+1))/2;
  ci = max(1, min(ndiv, round((vm + CLIM)/(2*CLIM)*(ndiv-1)) + 1));
  w = 1;
  if ~isempty(sg)
    sm = (sg(k) + sg(k+1))/2;
    if isfinite(sm) && sm > 0
      w = min(1, abs(vm)/(NSIG*sm));
      if abs(vm) >= NSIG*sm, nsig = nsig + 1; end
    end
  end
  col = mid + w*(dmap(ci,:) - mid);
  plot(ax, x(k:k+1), y(k:k+1), '-', 'Color', col, ...
    'LineWidth', lw, 'Clipping', 'on');
end
end
