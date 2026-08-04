function P = firnColumn(par)
%FIRNCOLUMN Density, refractive index and vertical twtt of the firn column.
%   P = FIRNCOLUMN(par) evaluates, on a depth grid d [m below the surface],
%
%     P.d     depth below the surface [m], ascending, column vector
%     P.rho   relative density rho/rho_ice, single-stage Herron-Langway
%     P.n     refractive index sqrt(eps_r) from the chosen mixing law
%     P.twtt  vertical (nadir) two-way traveltime from the surface [s]
%
%   The density profile is the single-stage Herron-Langway logistic
%
%     rho(d) = 1 / (1 + exp(beta + (log(1/rho_bco - 1) - beta)*d/bco_depth))
%
%   with beta = log(1/rho_sfc - 1), i.e. pinned exactly at rho_sfc at the
%   surface and rho_bco at bco_depth, approaching solid ice below. This is
%   the same functional form used by ptt.columnProfiles in the
%   fabric_anisotropy project, rewritten in depth-below-surface.
%
%   See also vdef.defaultParams, vdef.depthFromTwtt.

C = vdef.constants();

if nargin < 1 || isempty(par)
  par = vdef.defaultParams();
end

d = linspace(0, par.max_depth, par.nz).';

beta  = log(1/par.rho_sfc - 1);
slope = log(1/par.rho_bco - 1) - beta;   % change in logit over one bco_depth
rho   = 1 ./ (1 + exp(beta + slope * d / par.bco_depth));
rho   = min(rho, 1);

switch lower(par.mixing)
  case 'kovacs'
    % Kovacs et al. (1995), linear in density: n = 1 + a*rho[g/cm^3]
    n = 1 + C.kovacs_a * rho * (C.rho_ice/1000);
  case 'looyenga'
    % Looyenga (1965) / Glen-Paren: eps^(1/3) linear in relative density
    n = (1 + rho*(C.eps_ice^(1/3) - 1)).^(3/2);
  otherwise
    error('vdef:firnColumn:mixing', ...
      'Unknown mixing law "%s" (expected ''kovacs'' or ''looyenga'').', par.mixing);
end

% Vertical two-way traveltime below the surface (nadir, no ray bending)
P.d    = d;
P.rho  = rho;
P.n    = n;
P.twtt = (2/C.c) * cumtrapz(d, n);

end
