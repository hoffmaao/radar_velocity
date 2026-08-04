function [depth, n_local] = depthFromTwtt(P, twtt_below_surface)
%DEPTHFROMTWTT Map two-way traveltime below the surface to depth.
%   [depth, n_local] = DEPTHFROMTWTT(P, twtt) inverts the vertical twtt
%   table P (from vdef.firnColumn) for the depth below the surface and
%   returns the local refractive index at that depth.
%
%   twtt may be any shape; depth and n_local match it. Values outside the
%   table (negative twtt, or beyond P.max_depth) come back NaN, which
%   propagates as an invalid sample through the rest of the chain.
%
%   n_local is the quantity that converts a differential traveltime into a
%   physical displacement: a reflector whose two-way traveltime below the
%   surface changes by dtau has moved by dtau*c/(2*n_local) relative to the
%   surface, because the material added to or removed from the column is
%   the material AT that reflector's depth.
%
%   See also vdef.firnColumn, vdef.verticalDisplacement.

sz = size(twtt_below_surface);
t  = twtt_below_surface(:);

% P.twtt is strictly increasing (n > 0), so interp1 is well posed
depth   = interp1(P.twtt, P.d, t, 'linear', NaN);
n_local = interp1(P.d, P.n, depth, 'linear', NaN);

depth   = reshape(depth,   sz);
n_local = reshape(n_local, sz);

end
