function tide = cats2008_tide(gps_time)
%CATS2008_TIDE CATS2008 tide height at the survey site for given GPS times.
%   tide = CATS2008_TIDE(gps_time) interpolates the CATS2008 (v2023)
%   prediction in cats2008_apres_window.csv, which sits beside this file
%   and covers 7-13 Dec 2022 at 5-minute steps at the ApRES window on the
%   line, onto gps_time given in UTC seconds since 1970 (the OPR gps_time
%   convention). Times outside the window come back NaN.
%
%   WHY AN EXTERNAL TIDE. The surface admittance a(x) is a regression of
%   block heights on the tide across passes. Regressing on each pass's
%   own line-mean height instead - the chain's first estimator - puts
%   every per-pass height error into the regressor as well as the data,
%   and the slope is pulled toward one. See vdef.surfaceAdmittance for
%   the measured size of that error on the EAGER lines and what it does
%   to E*.
%
%   Tidal phase across the 5 km line is negligible at diurnal periods, so
%   one prediction serves every block. The prediction is used as a
%   REGRESSOR, not as an absolute datum: a scale or phase error in CATS
%   scales every block's admittance alike and cancels in the shape fit.

here = fileparts(mfilename('fullpath'));
fn = fullfile(here, 'cats2008_apres_window.csv');
assert(exist(fn, 'file') == 2, 'CATS2008 prediction not found: %s', fn);

% Plain csvread would choke on the ISO column; read the two numeric ones.
fid = fopen(fn, 'r');
C = textscan(fid, '%f %s %f', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
dn = C{1}; h = C{3};
assert(numel(dn) > 2 && all(isfinite(dn)) && all(isfinite(h)), ...
  'could not parse %s', fn);

t_dn = datenum(1970,1,1) + gps_time(:).' / 86400;
tide = interp1(dn, h, t_dn, 'linear', NaN);
end
