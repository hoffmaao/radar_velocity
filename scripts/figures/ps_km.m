function [xk, yk] = ps_km(ps, lon, lat)
%PS_KM Project lon/lat to EPSG:3031 Antarctic Polar Stereographic, in km.
[x, y] = projfwd(ps, lat, lon);
xk = x/1e3; yk = y/1e3;
end
