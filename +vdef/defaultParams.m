function par = defaultParams()
%DEFAULTPARAMS Default firn-ice column parameters for the twtt-depth map.
%
%   Coordinate convention: d is DEPTH BELOW THE SNOW SURFACE, in metres,
%   increasing downward. (Note this is the opposite of the +ptt package in
%   the fabric_anisotropy project, where zhat = z/H is height above the
%   bed. The density parameterisation is the same single-stage
%   Herron-Langway form, recast in depth and pinned at the surface and at
%   bubble close-off, so rho_sfc / rho_bco carry over unchanged.)
%
%   Fields:
%     rho_sfc    relative density (rho/rho_ice) of surface snow
%     rho_bco    relative density at bubble close-off
%     bco_depth  depth of bubble close-off below the surface [m]
%     mixing     'kovacs' (default) or 'looyenga': density -> refractive index
%     max_depth  extent of the depth grid [m]
%     nz         number of depth grid points
%
%   Only the density profile enters the vertical-deformation inversion, via
%   the refractive index n(d) that converts a differential two-way
%   traveltime into a physical displacement. At the accum3 band and the
%   depths of interest the result is insensitive to rho_sfc/rho_bco below
%   the firn-ice transition (n is pinned at n_ice there), so the defaults
%   only really matter for the top ~100 m.

par.rho_sfc   = 0.35;
par.rho_bco   = 0.81;
par.bco_depth = 60;

par.mixing    = 'kovacs';
par.max_depth = 1500;
par.nz        = 6001;

end
