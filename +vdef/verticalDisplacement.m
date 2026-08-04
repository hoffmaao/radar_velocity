function V = verticalDisplacement(blk, map, par, opts)
%VERTICALDISPLACEMENT Convert block-averaged dtau into vertical displacement.
%   V = VERTICALDISPLACEMENT(blk, map, par, opts) maps the block-averaged
%   differential traveltime onto physical depth and converts it into the
%   vertical displacement of each reflector RELATIVE TO THE SURFACE over
%   the repeat interval, and into the corresponding relative vertical
%   velocity.
%
%   THE CONVERSION. Let h(t) be the depth of a reflector below the surface.
%   Its vertical two-way traveltime is
%
%     tau(t) = (2/c) * int_0^{h(t)} n(z,t) dz
%
%   so between the two passes
%
%     dtau = (2/c) * [ n(h) * dh  +  int_0^h dn/dt dz * dt ]
%
%   The first term is the geometric one: the column between the surface and
%   the reflector has changed thickness by dh, and the material added or
%   removed is at the reflector's own depth, so the LOCAL index n(h) is the
%   right conversion factor - not a column average.
%
%     dh = dtau * c / (2*n(h))
%
%   The second term is the densification correction: firn compacting under
%   the reflector-free column raises n over time. Over the hours-to-days
%   baselines of the McMurdo repeat passes it is utterly negligible
%   (densification changes n by ~1e-5 per year in the upper firn), but it
%   grows linearly with the repeat interval and is not negligible for
%   season-to-season baselines. It is applied only when
%   opts.densification_rate is set to a non-zero value (relative density
%   per year, applied over the firn column above each reflector); the
%   default of 0 leaves it off and records that fact in V.densification_applied.
%
%   SIGN. dh < 0 means the column between the surface and the reflector has
%   SHORTENED over the interval, i.e. net compaction / vertical
%   compressive strain. dh > 0 means it has lengthened (vertical extension).
%
%   WHAT IS AND IS NOT MEASURED. Because dtau is referenced to the surface
%   return, everything common to the whole column drops out: antenna height
%   change, tidal heave of a floating column as a rigid body, bulk timing
%   drift, and the unwrapping constant. What survives is the DIFFERENTIAL
%   vertical motion within the ice column - which is exactly the vertical
%   strain signal - plus any tidal FLEXURE that strains the column, which
%   at a grounding-zone site like Windless Bight is a real and
%   time-varying part of the signal rather than an error.
%
%   Returns:
%     V.depth     Nt x Nblk depth below the surface [m]
%     V.n_local   Nt x Nblk local refractive index
%     V.dh        Nt x Nblk vertical displacement relative to surface [m]
%     V.dh_std    Nt x Nblk propagated 1-sigma of dh [m]. This is the error
%                 bar on V.dh: it comes from blk.dtau_std, the standard
%                 error of the block mean, not from the within-block
%                 scatter of the samples that formed it.
%     V.dh_scatter Nt x Nblk propagated within-block scatter of dh [m]. A
%                 relative noise measure for weighting, NOT an error bar.
%     V.v         Nt x Nblk relative vertical velocity [m/yr] (NaN if no dt)
%     V.v_std, V.v_scatter  the same two quantities as rates [m/yr]
%     V.delta_t   repeat interval used [yr]
%
%   See also vdef.blockAverage, vdef.invertStrainRate, vdef.firnColumn.

C = vdef.constants();
P = vdef.firnColumn(par);

Time = map.Time(:);
[Nt, Nblk] = size(blk.dtau);

% Depth of every (bin, block) from its twtt below that block's surface
twtt_below = bsxfun(@minus, Time, blk.Surface(:).');
[depth, n_local] = vdef.depthFromTwtt(P, twtt_below);

% Geometric conversion, using the LOCAL index at the reflector depth
dh     = blk.dtau     * C.c / 2 ./ n_local;
dh_std = blk.dtau_std * C.c / 2 ./ n_local;
if isfield(blk,'dtau_scatter')
  dh_scatter = blk.dtau_scatter * C.c / 2 ./ n_local;
else
  dh_scatter = dh_std;
end

% Optional densification correction (off by default)
densification_applied = false;
if isfield(opts,'densification_rate') && ~isempty(opts.densification_rate) ...
    && opts.densification_rate ~= 0 && isfield(opts,'delta_t') && ~isempty(opts.delta_t)
  % dn/dt = (dn/drho)*drho/dt, integrated from the surface to each depth.
  % Kovacs is linear in density so dn/drho is a constant; for looyenga the
  % local derivative is used.
  switch lower(par.mixing)
    case 'kovacs'
      dn_drho = C.kovacs_a * (C.rho_ice/1000) * ones(size(P.d));
    otherwise
      dn_drho = gradient(P.n, P.rho);
  end
  % Densification is confined to the firn: taper to zero at close-off
  taper = max(0, 1 - P.d/par.bco_depth);
  integrand = dn_drho .* opts.densification_rate .* taper;
  col = cumtrapz(P.d, integrand);                       % [1] per year, x depth
  corr_tau = (2/C.c) * interp1(P.d, col, depth, 'linear', NaN) * opts.delta_t;
  dh = dh - corr_tau * C.c / 2 ./ n_local;
  densification_applied = true;
end

V = [];
V.depth   = depth;
V.n_local = n_local;
V.dh      = dh;
V.dh_std  = dh_std;
V.dh_scatter = dh_scatter;
V.densification_applied = densification_applied;

if isfield(opts,'delta_t') && ~isempty(opts.delta_t) && opts.delta_t ~= 0
  V.delta_t   = opts.delta_t;
  V.v         = dh / opts.delta_t;
  V.v_std     = dh_std / abs(opts.delta_t);
  V.v_scatter = dh_scatter / abs(opts.delta_t);
else
  V.delta_t   = NaN;
  V.v         = nan(Nt, Nblk);
  V.v_std     = nan(Nt, Nblk);
  V.v_scatter = nan(Nt, Nblk);
end

end
