function drawn = overlay_gl(ax, gis_dir)
drawn = false;
fn = fullfile(gis_dir,'GroundingLine_Antarctica_v02.shp');
if ~exist(fn,'file'), return; end
try
  S = shaperead(fn);
  X = []; Y = [];
  for k = 1:numel(S)
    X = [X; S(k).X(:); NaN]; Y = [Y; S(k).Y(:); NaN]; %#ok<AGROW>
  end
  gx = X/1e3; gy = Y/1e3;      % the shapefile is already EPSG:3031 metres
  xl = xlim(ax); yl = ylim(ax); pad = 3;
  bad = gx < xl(1)-pad | gx > xl(2)+pad | gy < yl(1)-pad | gy > yl(2)+pad;
  gx(bad) = NaN; gy(bad) = NaN;
  if all(isnan(gx)), return; end
  plot(ax, gx, gy, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);
  xlim(ax, xl); ylim(ax, yl);
  drawn = true;
catch
end
end
