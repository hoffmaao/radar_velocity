function ok = rema_underlay(ax, tif)
%REMA_UNDERLAY Draw the REMA v2 hillshade under an EPSG:3031 axes in km.
%   The axes are already in the raster's own projection, so the block
%   covering the current view is read and dropped straight in - no
%   resampling, no rotation, no reprojection. It is lifted into a light
%   gray range so the data drawn on top stays dominant, and pushed to the
%   bottom of the draw order. The block is padded 8% beyond the limits so
%   a small later limit adjustment does not expose white strips. Returns
%   false (with a message) when the tile is missing or unreadable, and
%   the figure then renders exactly as it did without imagery.
%
%   The browse tif is a COG whose overview IFDs make geotiffinfo error
%   ("multiple images ... sizes are different"), so the georeference comes
%   from georasterinfo and the pixels from imread on IFD 1.
ok = false;
if ~exist(tif,'file')
  fprintf('REMA hillshade not found (%s) - no imagery underlay.\n', tif);
  return;
end
try
  R3 = georasterinfo(tif).RasterReference;
  xl = xlim(ax); yl = ylim(ax);                     % km, EPSG:3031
  xg = xl + 0.08*diff(xl)*[-1 1];
  yg = yl + 0.08*diff(yl)*[-1 1];
  px = R3.CellExtentInWorldX;
  py = R3.CellExtentInWorldY;                       % may differ from px
  c0 = max(1, floor((xg(1)*1e3 - R3.XWorldLimits(1))/px));
  c1 = min(R3.RasterSize(2), ceil((xg(2)*1e3 - R3.XWorldLimits(1))/px));
  r0 = max(1, floor((R3.YWorldLimits(2) - yg(2)*1e3)/py));
  r1 = min(R3.RasterSize(1), ceil((R3.YWorldLimits(2) - yg(1)*1e3)/py));
  if c1 <= c0 || r1 <= r0
    fprintf('REMA tile does not cover this view - no imagery underlay.\n');
    return;
  end
  A = double(imread(tif, 'Index', 1, 'PixelRegion', {[r0 r1],[c0 c1]}));
  A(A == 0) = NaN;                                  % nodata
  vv = sort(A(isfinite(A)));
  if isempty(vv)
    fprintf('REMA tile is empty over this view - no imagery underlay.\n');
    return;
  end
  lo = vv(max(1,round(0.02*numel(vv)))); hi = vv(round(0.98*numel(vv)));
  g = min(max((A - lo)/max(hi-lo, 1), 0), 1);
  g = 0.58 + 0.40*g;                                % recessive light grays
  g(~isfinite(A)) = 1;                              % nodata as paper white
  % cell centres of the block actually read, in km; y descends with row
  ximg = (R3.XWorldLimits(1) + ([c0 c1] - 0.5)*px)/1e3;
  yimg = (R3.YWorldLimits(2) - ([r0 r1] - 0.5)*py)/1e3;
  hImg = image(ax, 'XData', ximg, 'YData', yimg, 'CData', repmat(g,[1 1 3]));
  uistack(hImg, 'bottom');
  xlim(ax, xl); ylim(ax, yl);
  ok = true;
  fprintf('REMA hillshade underlay: %d x %d px, native EPSG:3031\n', ...
    r1-r0+1, c1-c0+1);
catch ME
  fprintf('REMA underlay failed (%s) - continuing without imagery.\n', ME.message);
end
end
