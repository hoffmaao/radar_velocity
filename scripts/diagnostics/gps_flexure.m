%GPS_FLEXURE Tidal flexure profile from the GPS platform heights alone.
%
% WHY. The radar admittance profile correlates at -0.84 to -0.97 with an
% artefact predictor built from ref_z. That is not proof the hinge is
% spurious, because ref_z is confounded: on a floating shelf the platform
% rides the ice, so the differential platform height between two passes IS
% the differential tidal displacement, which carries the REAL flexure. The
% artefact and the signal are driven by the same physical quantity.
%
% This measures the flexure directly from the GPS, with no radar at all.
% For along-track position x and pass k, ref_z_k(x) is the platform height
% relative to the main pass, already on the main pass's along-track axis.
% Take each pass's LINE MEAN as its tide, then per position regress
% ref_z_k(x) on that tide across passes:
%
%   ref_z_k(x) = a(x) * tide_k + const
%
% a(x) is the local vertical tidal admittance of the ice surface: ~1 where
% the shelf floats freely, falling toward 0 into the grounding zone. The
% hinge is where a(x) drops. If the GPS puts the hinge where the radar
% does, real flexure IS present there whatever the radar contamination.

mp_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
BLOCK  = 200;
names  = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2','EAGER_2022_GL3','EAGER_2022_GL4'};

for n = 1:numel(names)
  pn = names{n};
  L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass','param_multipass');
  main_idx = L.param_multipass.multipass.baseline_master_idx;
  Np = numel(L.pass);
  Nx = numel(L.pass(main_idx).ref_z);

  Z = nan(Nx, Np); tide = nan(1,Np);
  for k = 1:Np
    z = L.pass(k).ref_z(:);
    if numel(z) ~= Nx, continue; end
    Z(:,k) = z;
    tide(k) = mean(z,'omitnan');
  end
  lat = L.pass(main_idx).lat(:); lon = L.pass(main_idx).lon(:);
  alongt = L.pass(main_idx).along_track(:);
  clear L;

  nb = floor(Nx/BLOCK);
  a = nan(nb,1); xb = nan(nb,1); latb = nan(nb,1);
  for b = 1:nb
    idx = (b-1)*BLOCK+1 : b*BLOCK;
    zb = mean(Z(idx,:), 1, 'omitnan');
    ok = isfinite(zb) & isfinite(tide);
    if nnz(ok) < 5, continue; end
    p = polyfit(tide(ok), zb(ok), 1);
    a(b) = p(1);
    xb(b) = mean(alongt(idx),'omitnan');
    latb(b) = mean(lat(idx),'omitnan');
  end

  fprintf('\n===== %s (main %d, %d passes) =====\n', pn, main_idx, Np);
  fprintf('tide range across passes: %.3f m\n', max(tide)-min(tide));
  fprintf('%6s %10s %10s %9s\n','block','along[km]','lat','a(x)');
  for b = 1:nb
    if ~isfinite(a(b)), continue; end
    fprintf('%6d %10.2f %10.4f %9.3f\n', b, xb(b)/1e3, latb(b), a(b));
  end
  % change point in a(x), same detector the radar analysis uses
  ok = find(isfinite(a)); MIN_SIDE = 3;
  if numel(ok) >= 2*MIN_SIDE
    av = a(ok); best = -inf; kb = [];
    for k = MIN_SIDE:(numel(ok)-MIN_SIDE)
      dd = mean(av(k+1:end)) - mean(av(1:k));   % a should RISE seaward
      if dd > best, best = dd; kb = k; end
    end
    fprintf('largest rise in a(x): %.3f at along %.2f km (lat %.4f)\n', ...
      best, 0.5*(xb(ok(kb))+xb(ok(kb+1)))/1e3, 0.5*(latb(ok(kb))+latb(ok(kb+1))));
  end
end
