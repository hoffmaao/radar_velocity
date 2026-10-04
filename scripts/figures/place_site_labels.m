function place_site_labels(ax, AP, ink)
%PLACE_SITE_LABELS Site names beside their markers, pushed clear of each other.
%   GA04 and GA05 are ~0.1 km apart, so labels anchored at the marker
%   overlap and neither reads. Sites are taken in y order and each label is
%   pushed up until it clears the previous one; a label moved far enough
%   that its pairing stops being obvious gets a hairline leader back to its
%   own marker. Call AFTER the limits are final - the separation and the
%   offset are both fractions of the axis extent.
if isempty(AP), return; end
xl = xlim(ax); yl = ylim(ax);
sep = 0.030*diff(yl);          % minimum label separation, data units
dx  = 0.045*diff(xl);          % label offset to the right of the marker
[~, ord] = sort([AP.y]);
ylab = nan(1, numel(AP)); prev = -inf;
for q = ord(:).'
  ylab(q) = max(AP(q).y/1e3, prev + sep);
  prev = ylab(q);
end
for q = 1:numel(AP)
  xq = AP(q).x/1e3; yq = AP(q).y/1e3;
  if abs(ylab(q) - yq) > 0.3*sep
    plot(ax, [xq + 0.25*dx, xq + 0.85*dx], [yq, ylab(q)], '-', ...
      'Color', ink, 'LineWidth', 0.5);
  end
  text(ax, xq + dx, ylab(q), AP(q).name, 'FontSize', 8, 'Color', ink, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','left');
end
end
