%BETWEEN_BUILD Is the 8.28 mm floor a property of the method, or of one product?
%
%   The error budget's dominant term - the 8.28 mm between-build
%   disagreement - was measured on the one comparison where two products
%   image the SAME ice: EAGER_2022 vs GL1. But EAGER_2022 is also the one
%   wholly uncalibrated product (coregistration_time_shift and
%   equalization both absent), so that number may characterise one bad
%   build rather than the method.
%
%   This compares EVERY product pair. The legs are 120-145 m apart, which
%   against a >=1 km flexure scale means the true tidal admittance should
%   differ by little between adjacent legs - so cross-leg disagreement is
%   an upper bound on the method's reproducibility that does not involve
%   EAGER_2022. The matrix separates the hypotheses:
%     - every pair disagrees at ~8 mm     -> the floor is the method
%     - only EAGER_2022 rows are ~8 mm    -> the floor is that product,
%       and the real reproducibility is whatever GL1..GL4 show
%
%   Both quantities are reported per pair: tidal admittance [mm per m of
%   tide] and secular change over the window [mm], top 100 m, plus the
%   expected disagreement from the two fits' own sigmas, so excess over
%   formal error is visible directly.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../between_build.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

if ~exist('VVEL_SUFFIX','var'), VVEL_SUFFIX = '_v3'; end
vvel_dir = ['/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel' VVEL_SUFFIX];
mp_dir   = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';

PASS_NAMES = {'EAGER_2022','EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
REF_DEPTH    = 100;
MAX_BASELINE = 10;
TOL          = 100;    % [m] blocks closer than this along track are compared

%% Fit every product
R = [];
for n = 1:numel(PASS_NAMES)
  S = one_product(PASS_NAMES{n}, vvel_dir, mp_dir, REF_DEPTH, MAX_BASELINE);
  if isempty(S), continue; end
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
end
assert(numel(R) >= 2, 'need at least two products');

%% Pairwise disagreement matrices
np = numel(R);
fprintf('\n===== between-build disagreement, top %d m =====\n', REF_DEPTH);
fprintf('%-18s %-18s %7s | %10s %10s %8s | %10s %10s\n', ...
  'product A','product B','xsep[m]','adm rms','adm expct','excess', ...
  'sec rms','sec expct');
worst_e = 0; worst_c = 0;
for a = 1:np
  for b = a+1:np
    [adm_rms, adm_exp, sec_rms, sec_exp, nblk] = compare_pair(R(a), R(b), TOL);
    if ~isfinite(adm_rms), continue; end
    xsep = abs(R(a).xtrack - R(b).xtrack);
    exc = adm_rms/max(adm_exp,eps);
    fprintf('%-18s %-18s %7.0f | %10.2f %10.2f %7.1fx | %10.2f %10.2f\n', ...
      R(a).name, R(b).name, xsep, adm_rms, adm_exp, exc, sec_rms, sec_exp);
    inv_e = any(strcmp({R(a).name,R(b).name},'EAGER_2022'));
    if inv_e, worst_e = max(worst_e, adm_rms); else, worst_c = max(worst_c, adm_rms); end
  end
end
fprintf('\nworst admittance rms involving EAGER_2022:  %.2f mm\n', worst_e);
fprintf('worst admittance rms among calibrated only: %.2f mm\n', worst_c);
fprintf(['\nReading: if the calibrated-only worst case sits far below the\n' ...
  'EAGER_2022 rows, the 8.28 mm floor was that product, not the method,\n' ...
  'and the calibrated number is the honest reproducibility floor.\n']);

%% ========================================================================
function S = one_product(pn, vvel_dir, mp_dir, REF_DEPTH, MAX_BASELINE)
S = [];
f = dir(fullfile(vvel_dir, [pn '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f), warning('no products for %s', pn); return; end

L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
Np = numel(L.pass); elev = nan(1,Np);
xtr = nan(1,Np);
for k = 1:Np
  elev(k) = mean(L.pass(k).elev,'omitnan');
  if isfield(L.pass(k),'ref_y'), xtr(k) = median(L.pass(k).ref_y,'omitnan'); end
end
clear L;

spy = 365.25*86400; np = numel(f);
strain = []; tide = nan(1,np); tday = nan(1,np); mb = nan(1,np);
along = []; ref0 = NaN;
for i = 1:np
  o = load(fullfile(vvel_dir, f(i).name));
  if isnan(ref0), ref0 = o.pass_idx_ref; end
  if isempty(strain)
    strain = nan(numel(o.S1), np);
    along = o.Along_track(:);
  end
  for b = 1:size(strain,1)
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    strain(b,i) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
  end
  tide(i) = elev(o.pass_idx_sec) - elev(ref0);
  tday(i) = mean(o.GPS_time + o.delta_t_blk*spy,'omitnan')/86400;
  mb(i)   = max(abs(o.baseline_y));
end
use = ~(isfinite(mb) & mb > MAX_BASELINE);
if nnz(use) < 6, return; end
strain = strain(:,use); tide = tide(use); tday = tday(use) - min(tday(use));
span = max(tday) - min(tday);

A = vdef.fitTideAdmittance(strain, tday, tide);
S = struct('name',pn,'along',along, ...
  'adm', 1e3*REF_DEPTH*A.admittance(:), ...
  'adm_std', 1e3*REF_DEPTH*A.admittance_std(:), ...
  'sec', 1e3*REF_DEPTH*A.trend(:)*span, ...
  'sec_std', 1e3*REF_DEPTH*A.trend_std(:)*span, ...
  'xtrack', median(xtr,'omitnan'));
end

%% ========================================================================
function [adm_rms, adm_exp, sec_rms, sec_exp, n] = compare_pair(A, B, TOL)
adm_rms = NaN; adm_exp = NaN; sec_rms = NaN; sec_exp = NaN; n = 0;
da = []; ea = []; ds = []; es = [];
for k = 1:numel(A.along)
  [gap, m] = min(abs(B.along - A.along(k)));
  if gap > TOL, continue; end
  if isfinite(A.adm(k)) && isfinite(B.adm(m))
    da(end+1) = A.adm(k) - B.adm(m); %#ok<AGROW>
    ea(end+1) = A.adm_std(k)^2 + B.adm_std(m)^2; %#ok<AGROW>
  end
  if isfinite(A.sec(k)) && isfinite(B.sec(m))
    ds(end+1) = A.sec(k) - B.sec(m); %#ok<AGROW>
    es(end+1) = A.sec_std(k)^2 + B.sec_std(m)^2; %#ok<AGROW>
  end
end
n = numel(da);
if n < 4, return; end
adm_rms = sqrt(mean(da.^2)); adm_exp = sqrt(mean(ea));
sec_rms = sqrt(mean(ds.^2)); sec_exp = sqrt(mean(es));
end
