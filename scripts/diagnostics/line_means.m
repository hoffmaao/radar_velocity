%LINE_MEANS Line-mean tidal response per product, and does block size matter?
%
%   TWO QUESTIONS, one script.
%
%   1. ALONG-TRACK AVERAGING. The block is the averaging unit: 200 columns
%      = 500 m at standard processing. Averaging harder (1000 columns =
%      2.5 km) buys more looks per block at the cost of along-track
%      resolution - which for a LINE MEAN we do not need. If the per-block
%      noise is incoherent speckle, the line mean should come out the same
%      with a similar sigma either way (same data, rearranged); if some
%      per-block error is nonlinear (unwrap, coverage edges), coarse
%      blocks can genuinely beat fine ones. Measured, not argued.
%
%   2. THE CROSS-PRODUCT MEAN. The between-build matrix showed calibrated
%      builds agree at 1.2-1.5x their formal errors, which makes averaging
%      their line means defensible. The last table pools GL1-GL4
%      (EAGER_2022 is excluded twice over: it is a second build of GL1
%      rather than an independent line - see vdef.surveyLines - and it is
%      wholly uncalibrated, 2.6x excess) with
%      inverse-variance weights, inflated by the measured 1.4x
%      between-build excess so the quoted sigma is honest rather than
%      formal.
%
%   Everything from the network inversion over all pairs - no reference
%   pass anywhere.
%
%   Run on the server after chain2 (needs CSARP_vvel_net and _net1k):
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../line_means.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root   = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
REF_DEPTH = 100; MAX_BASELINE = 10;
EXCESS = 1.4;         % measured between-build excess factor, between_build.m
SIG_APRES = -1.24;    % mm per m of tide over the top 100 m, ApRES GA04
                      % rate method (+/- 0.04), apres_rate_check.py

% The four lines, plus the duplicate build of GL1. EAGER_2022 is NOT a
% fifth profile - it is a second build of GL1 sharing thirteen of its
% fourteen passes (vdef.surveyLines) - but the floor computed below needs
% two builds of the SAME ice, so it is loaded here deliberately and must
% not be read as another line.
[LINES, DUP] = vdef.surveyLines();
PASS_NAMES = [{DUP.name}, LINES];
DIRS = {fullfile(root,'CSARP_vvel_net'), '500 m blocks'; ...
        fullfile(root,'CSARP_vvel_net1k'), '2.5 km blocks'};

res = nan(numel(PASS_NAMES), size(DIRS,1), 4);   % [adm, adm_sig, sec, sec_sig]
for di = 1:size(DIRS,1)
  fprintf('\n===== %s (%s) =====\n', DIRS{di,2}, DIRS{di,1});
  fprintf('%-18s %6s %7s | %10s %9s | %10s %9s\n','line','pairs','blocks', ...
    'adm[mm/m]','1-sigma','sec[mm]','1-sigma');
  for n = 1:numel(PASS_NAMES)
    [adm, adm_sig, sec, sec_sig, npair, nblk] = ...
      line_mean(PASS_NAMES{n}, DIRS{di,1}, mp_dir, REF_DEPTH, MAX_BASELINE);
    if ~isfinite(adm)
      fprintf('%-18s SKIP (products missing or too few pairs)\n', PASS_NAMES{n});
      continue;
    end
    res(n,di,:) = [adm adm_sig sec sec_sig];
    fprintf('%-18s %6d %7d | %+10.2f %9.2f | %+10.2f %9.2f\n', ...
      PASS_NAMES{n}, npair, nblk, adm, adm_sig, sec, sec_sig);
  end
end

%% Block-size verdict
fprintf('\n===== 500 m vs 2.5 km blocks, line means =====\n');
fprintf('%-18s %12s %12s %10s\n','line','adm@500m','adm@2.5km','diff/sig');
for n = 1:numel(PASS_NAMES)
  a = res(n,1,:); b = res(n,2,:);
  if ~isfinite(a(1)) || ~isfinite(b(1)), continue; end
  z = (a(1)-b(1))/sqrt(a(2)^2 + b(2)^2);
  fprintf('%-18s %+9.2f/%.2f %+9.2f/%.2f %9.1f\n', PASS_NAMES{n}, ...
    a(1), a(2), b(1), b(2), z);
end

%% Cross-product mean, calibrated lines only
fprintf('\n===== pooled line mean, GL1-GL4, both block sizes =====\n');
for di = 1:size(DIRS,1)
  v = squeeze(res(2:5,di,1)); s = squeeze(res(2:5,di,2)) * EXCESS;
  ok = isfinite(v) & isfinite(s) & s > 0;
  if nnz(ok) < 2, continue; end
  w = 1./s(ok).^2;
  mu = sum(w.*v(ok))/sum(w); se = sqrt(1/sum(w));
  % scatter-based check on the quoted error
  se_scat = std(v(ok))/sqrt(nnz(ok));
  fprintf('%-14s: %+6.2f +/- %.2f mm per m of tide (scatter-based %.2f; ApRES rate method %+.2f, %.1f sigma apart)\n', ...
    DIRS{di,2}, mu, se, se_scat, SIG_APRES, abs(mu - SIG_APRES)/max(se,se_scat));
end
fprintf(['\nThe pooled value is the project''s best radar estimate of the tidal\n' ...
  'response. Sigma includes the measured %.1fx between-build excess.\n'], EXCESS);

%% ========================================================================
function [adm, adm_sig, sec, sec_sig, npair, nblk] = ...
    line_mean(pn, vdir, mdir, REF_DEPTH, MAX_BASELINE)
adm = NaN; adm_sig = NaN; sec = NaN; sec_sig = NaN; npair = 0; nblk = 0;
L = load(fullfile(mdir, sprintf('%s_multipass03.mat', pn)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
clear L;

f = dir(fullfile(vdir, [pn '_vvel_*.mat']));
P = []; D = []; W = [];
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  if isempty(tok), continue; end
  o = load(fullfile(vdir, f(q).name));
  if ~vdef.pairAligned(o), continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  Nblk = numel(o.S1); sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
  end
  if all(~isfinite(sv)), continue; end
  if isempty(D), D = sv; else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < 10, return; end
npair = size(P,1); nblk = size(D,1);

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));
tday = (ptime - min(ptime))/86400;
tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);
span = max(tday) - min(tday);

va = 1e3*REF_DEPTH*A.admittance(:); sa = 1e3*REF_DEPTH*A.admittance_std(:);
vs = 1e3*REF_DEPTH*A.trend(:)*span; ss = 1e3*REF_DEPTH*A.trend_std(:)*span;
ok = isfinite(va) & isfinite(sa) & sa > 0;
if nnz(ok) < 2, return; end
w = 1./sa(ok).^2;
adm = sum(w.*va(ok))/sum(w);
% quoted sigma: the larger of formal and block-scatter, so few coarse
% blocks cannot fake precision through a lucky agreement
adm_sig = max(sqrt(1/sum(w)), std(va(ok))/sqrt(nnz(ok)));
ok = isfinite(vs) & isfinite(ss) & ss > 0;
w = 1./ss(ok).^2;
sec = sum(w.*vs(ok))/sum(w);
sec_sig = max(sqrt(1/sum(w)), std(vs(ok))/sqrt(nnz(ok)));
end
