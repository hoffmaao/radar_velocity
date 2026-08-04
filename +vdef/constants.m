function C = constants()
%CONSTANTS Physical constants for the vertical-deformation chain.
%
%   Units are SI throughout the +vdef package (metres, seconds, Hz), with
%   the single exception that reported strain rates and velocities are
%   converted to per-year / metres-per-year at the reporting stage only.
%
%   The Kovacs coefficient is chosen so that solid ice reproduces the OPR
%   physical_constants value er_ice = 3.15:
%     n_ice = 1 + 0.845*0.917 = 1.7749,  n_ice^2 = 3.150

C.c        = 299792458;    % speed of light in vacuum [m/s]
C.rho_ice  = 917;          % density of solid ice [kg/m^3]
C.eps_ice  = 3.15;         % relative permittivity of solid ice (OPR er_ice)
C.kovacs_a = 0.845;        % Kovacs et al. (1995): n = 1 + a*rho [rho in g/cm^3]
C.sec_per_year = 365.25*86400;

end
