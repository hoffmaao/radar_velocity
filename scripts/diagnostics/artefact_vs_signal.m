%ARTEFACT_VS_SIGNAL Is the tidal response explained by the residual artefact?
%
%   THIS IS THE ACCEPTANCE TEST FOR ANY TIDAL RESULT. It is what showed
%   that the flexure hinge reported before 2026-08-05 was the residual
%   misalignment and not ice. Run it before believing any along-track
%   structure in the admittance.
%
%   Results to date, corr(admittance, artefact predictor) per product:
%     scalar coalignment (CSARP_vvel_v2):
%       -0.97 / -0.93 / -0.84 / -0.95 / -0.92, every one p < 0.005
%     per-column coalignment (CSARP_vvel_v3):
%       +0.34 / -0.63 / +0.29 / +0.55 / -0.64, none significant
%   and the hinge disappeared from all five products at the same time.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../artefact_vs_signal.m')"

%
% Both coregistration_time_shift and coalignPair remove a SCALAR per pair.
% What no scalar can remove is the ALONG-TRACK VARIATION of the
% misalignment, so the scalar-invariant measure of what survives is
%
%   delta(x) = alpha * -(ref_z_sec(x) - ref_z_ref(x))/(c/2)
%   resid(x) = delta(x) - mean(delta)          [line mean removed]
%
% averaged into the same along-track blocks the vvel run used. (An earlier
% version of this diagnostic subtracted dtau_bulk instead of the line mean,
% which double-counted coregistration_time_shift on the coregistered
% products and made them look worse than the uncoregistered ones.)
%
% THE TEST. Per block, regress resid against tide across pairs: the slope
% g(b) [ns per metre of tide] is how much tide-proportional misalignment
% that block still carries - a per-block ARTEFACT PREDICTOR. If the tidal
% admittance profile is really the artefact, admittance(b) must track
% g(b) across blocks. If the hinge is ice, it must not.

mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
if ~exist('VVEL_SUFFIX','var'), VVEL_SUFFIX = '_v3'; end
vvel_dir = ['/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel' VVEL_SUFFIX];
addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef
c = 299792458; ALPHA = 1.03; BLOCK = 200; REF_DEPTH = 100; MAX_BASELINE = 10;

names = vdef.surveyLines();   % four lines; EAGER_2022 duplicates GL1

fprintf('%-16s %7s %11s %11s %13s %11s\n', 'product','pairs', ...
  'medRMS[ns]','max|g|[ns/m]','corr(adm,g)','p');
for n = 1:numel(names)
  pn = names{n};
  f = dir(fullfile(vvel_dir, [pn '_vvel_*.mat']));
  keep = ~cellfun('isempty', regexp({f.name}, ...
    ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
  f = f(keep);
  if isempty(f), continue; end

  L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
  Np = numel(L.pass); elev = nan(1,Np);
  for k = 1:Np, elev(k) = mean(L.pass(k).elev,'omitnan'); end

  sec_per_year = 365.25*86400;
  resid = []; strain = []; tide = nan(1,numel(f)); tday = nan(1,numel(f));
  maxbl = nan(1,numel(f));
  for i = 1:numel(f)
    o = load(fullfile(vvel_dir, f(i).name));
    zs = L.pass(o.pass_idx_sec).ref_z(:);
    zr = L.pass(o.pass_idx_ref).ref_z(:);
    d  = ALPHA * -(zs - zr)/(c/2);
    d  = d - mean(d,'omitnan');            % scalar-invariant: line mean out
    Nx = numel(d); nb = floor(Nx/BLOCK);
    if isempty(resid)
      Nblk = numel(o.S1);
      resid  = nan(Nblk, numel(f));
      strain = nan(Nblk, numel(f));
    end
    for b = 1:min(nb, size(resid,1))
      resid(b,i) = mean(d((b-1)*BLOCK+1 : b*BLOCK), 'omitnan');
    end
    for b = 1:size(strain,1)
      dd = o.depth_blk(:,b); ok = isfinite(dd) & isfinite(o.dh_blk(:,b));
      if ~any(ok) || max(dd(ok)) < REF_DEPTH, continue; end
      strain(b,i) = interp1(dd(ok), o.dh_blk(ok,b), REF_DEPTH, 'linear', NaN)/REF_DEPTH;
    end
    tide(i)  = elev(o.pass_idx_sec) - elev(o.pass_idx_ref);
    tday(i)  = mean(o.GPS_time + o.delta_t_blk*sec_per_year,'omitnan')/86400;
    maxbl(i) = max(abs(o.baseline_y));
  end
  clear L;

  use = ~(isfinite(maxbl) & maxbl > MAX_BASELINE);
  resid = resid(:,use); strain = strain(:,use);
  tide = tide(use); tday = tday(use) - min(tday(use));

  % per-block artefact predictor: ns of residual per metre of tide
  Nblk = size(resid,1); g = nan(Nblk,1);
  for b = 1:Nblk
    v = resid(b,:); ok = isfinite(v) & isfinite(tide);
    if nnz(ok) < 5, continue; end
    p = polyfit(tide(ok), v(ok), 1); g(b) = p(1)*1e9;
  end

  A = vdef.fitTideAdmittance(strain, tday, tide);
  adm = A.admittance(:);

  ok = isfinite(adm) & isfinite(g);
  rr = NaN; pv = NaN;
  if nnz(ok) >= 5
    m = corrcoef(g(ok), adm(ok)); rr = m(1,2);
    nn = nnz(ok); tstat = rr*sqrt((nn-2)/max(1-rr^2,eps));
    pv = betainc((nn-2)/((nn-2)+tstat^2), (nn-2)/2, 0.5);
  end
  fprintf('%-16s %7d %11.2f %12.2f %13.2f %11.3f\n', pn, nnz(use), ...
    median(sqrt(mean(resid.^2,2,'omitnan'))*1e9,'omitnan'), ...
    max(abs(g)), rr, pv);
end

fprintf('\ncorr(adm,g): correlation across blocks between the measured tidal\n');
fprintf('admittance and the residual-misalignment artefact predictor.\n');
fprintf('Near zero => the along-track profile is NOT the artefact.\n');
