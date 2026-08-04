function F = forwardDisplacement(eps_zz_fn, twtt, Surface, par, delta_t)
%FORWARDDISPLACEMENT Synthetic dtau for a prescribed vertical strain rate.
%   F = FORWARDDISPLACEMENT(eps_zz_fn, twtt, Surface, par, delta_t) is the
%   forward counterpart of the inversion chain: given a vertical strain
%   rate profile it produces the differential two-way traveltime a
%   repeat-pass interferometer would measure.
%
%     eps_zz_fn  function handle eps_zz(d) [1/yr], d = depth below surface
%     twtt       Nt x 1 fast-time axis [s]
%     Surface    scalar or 1 x Nx surface twtt [s]
%     par        firn column parameters (vdef.defaultParams)
%     delta_t    repeat interval [yr]
%
%   The chain inverted by the package is
%
%     eps_zz(d) -> v(d) = int_0^d eps_zz  -> dh = v*delta_t
%               -> dtau = 2*n(d)*dh/c
%
%   and this function walks it forward exactly, so a round trip through
%   vdef.differentialRange / blockAverage / verticalDisplacement /
%   invertStrainRate recovers the input to within the noise added by the
%   caller. It is the basis of scripts/synthetic_vertical_velocity.m and of
%   opr_vvel/test/test_vvel_task.m.
%
%   Returns:
%     F.dtau   Nt x Nx differential traveltime [s], NaN above the surface
%     F.depth  Nt x Nx depth below the surface [m]
%     F.v      Nt x Nx prescribed relative vertical velocity [m/yr]
%     F.dh     Nt x Nx prescribed displacement [m]
%
%   See also vdef.verticalDisplacement, vdef.invertStrainRate.

C = vdef.constants();
P = vdef.firnColumn(par);

twtt = twtt(:);
Nt = numel(twtt);
Surface = Surface(:).';
Nx = numel(Surface);

% v(d) = integral of the strain rate from the surface down to d
v_tab = cumtrapz(P.d, eps_zz_fn(P.d));

twtt_below = bsxfun(@minus, twtt, Surface);
[depth, n_local] = vdef.depthFromTwtt(P, twtt_below);

v  = interp1(P.d, v_tab, depth, 'linear', NaN);
dh = v * delta_t;

dtau = 2 * n_local .* dh / C.c;
dtau(twtt_below < 0) = NaN;

F = [];
F.dtau  = dtau;
F.depth = depth;
F.v     = v;
F.dh    = dh;

end
