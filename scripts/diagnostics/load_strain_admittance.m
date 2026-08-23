function S = load_strain_admittance(pn, net_dir, mp_dir, opts)
%LOAD_STRAIN_ADMITTANCE Per-block tidal admittance of englacial strain.
%   S = LOAD_STRAIN_ADMITTANCE(pn, net_dir, mp_dir, opts) loads every
%   all-pairs vvel product for line pn, inverts the pair network for
%   per-pass column change, and fits the reference-invariant joint model
%   dh = a + b*t + c*tide per along-track block. It is the same chain as
%   scripts/figures/admittance_map.m uses - same pair filters, same
%   network inversion, same admittance fit - factored out so the flexure
%   driver and the map cannot drift onto different estimators.
%
%   The returned admittance is in SI: METRES of column change at
%   opts.ref_depth per METRE of tide (the map script scales to mm at the
%   last moment; this one leaves scaling to the caller).
%
%   opts fields, all optional:
%     .ref_depth     depth the column change is read at [m] (default 100)
%     .max_baseline  reject pairs with |cross-track baseline| above this
%                    [m] (default 10)
%     .min_pairs     fewest usable pairs to attempt the line (default 10)
%
%   Returns [] when the line has no usable products, else:
%     S.along    block centres on the main-pass along-track axis [m]
%     S.lat, S.lon
%     S.adm      dh at ref_depth per metre of tide [m/m]
%     S.adm_std  1-sigma from the joint fit [m/m]
%     S.n_pair   pairs used
%
%   See also scripts/figures/admittance_map.m, vdef.invertNetwork,
%   vdef.fitTideAdmittance.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'ref_depth')    || isempty(opts.ref_depth),    opts.ref_depth = 100;   end
if ~isfield(opts,'max_baseline') || isempty(opts.max_baseline), opts.max_baseline = 10; end
if ~isfield(opts,'min_pairs')    || isempty(opts.min_pairs),    opts.min_pairs = 10;    end

S = [];
f = dir(fullfile(net_dir, [pn '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f)
  fprintf('  %-5s no vvel products under %s\n', strrep(pn,'EAGER_2022_',''), net_dir);
  return;
end

L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
clear L;

P = []; D = []; W = []; along = []; lat = []; lon = []; Nblk = 0;
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  o = load(fullfile(net_dir, f(q).name));
  if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
  if max(abs(o.baseline_y)) > opts.max_baseline, continue; end
  if Nblk == 0
    Nblk = numel(o.S1); along = o.Along_track(:);
    lat = o.Latitude(:); lon = o.Longitude(:);
  end
  sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < opts.ref_depth, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), opts.ref_depth, 'linear', NaN);
  end
  if all(~isfinite(sv)), continue; end
  if isempty(D), D = sv; else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < opts.min_pairs
  fprintf('  %-5s only %d usable pairs - strain admittance not built\n', ...
    strrep(pn,'EAGER_2022_',''), size(P,1));
  return;
end

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));
tday = (ptime - min(ptime))/86400; tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);

S = struct('along', along, 'lat', lat, 'lon', lon, ...
  'adm', A.admittance(:), 'adm_std', A.admittance_std(:), ...
  'n_pair', size(P,1), 'ref_depth', opts.ref_depth);
end
