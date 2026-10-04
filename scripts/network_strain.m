%NETWORK_STRAIN Strain from the FULL pair network instead of one reference.
%
%   THE POINT. The 'main' pairing measures every pass against one chosen
%   reference: 12 pairs out of the 78 that 13 passes admit. That discards
%   85% of the available measurements, makes the answer depend on which
%   pass was picked as reference, and leaves no internal check. Loop
%   closure showed the cost - the direct and chained paths disagreed by
%   0.96 to 3.23 mm rms, as large as the signal.
%
%   This inverts the whole network instead (vdef.invertNetwork): solve
%   d(p) = x(j) - x(i) in least squares for the per-pass displacement x,
%   with a sum(x) = 0 datum so no epoch is privileged, and with robust
%   rejection so inconsistent pairs identify themselves through the fit
%   residuals rather than through a proxy quality metric. Then fit
%   x = a + b*t + c*tide exactly as before.
%
%   Three things should improve, and the script reports all three:
%     - precision, because 78 measurements constrain 13 unknowns
%     - reference-independence, by construction
%     - honesty, because the residuals are the closure errors
%
%   Requires a product processed with pairs='all':
%     matlab -batch "only_pass_names={'EAGER_2022_GL3'}; pairing_override='all'; out_suffix='_net'; run('.../run_vvel_scratch.m')"
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../network_strain.m')"

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net');
ref_dir = fullfile(root,'CSARP_vvel_v3');
if ~exist('PASS_NAME','var'), PASS_NAME = 'EAGER_2022_GL3'; end
REF_DEPTH    = 100;
MAX_BASELINE = 10;

%% Pass metadata
L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', PASS_NAME)), 'pass');
Np = numel(L.pass);
elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end
clear L;

%% Load every pair of the network
[P, D, W, along, nskip] = load_network(net_dir, PASS_NAME, REF_DEPTH, MAX_BASELINE);
assert(~isempty(P), 'no network pairs found in %s', net_dir);
Nblk = size(D,1);
fprintf('%s: %d pairs loaded, %d skipped (baseline cut or coalignment not applied)\n', ...
  PASS_NAME, size(P,1), nskip);
fprintf('  %d passes, %d blocks; a single-reference pairing would use %d pairs\n', ...
  Np, Nblk, Np-1);

%% Invert
N = vdef.invertNetwork(P, D, struct('n_sigma', 3, 'weights', W, 'n_epoch', Np));
fprintf('  network: median %d of %d pairs kept, %.1f epochs solved per block\n', ...
  round(median(N.n_used)), size(P,1), mean(N.n_epoch));
fprintf('  closure (fit residual) rms: median %.3f mm, worst block %.3f mm\n', ...
  1e3*REF_DEPTH*median(N.rms,'omitnan'), 1e3*REF_DEPTH*max(N.rms));

% pairs rejected most often - the network naming its own bad data. Only
% rows the inversion actually solved can reject a pair; unsolved rows
% leave N.used all-false without meaning rejection.
solved = any(isfinite(N.x), 2);
rej = sum(~N.used(solved,:) & isfinite(D(solved,:)), 1);
[~, ord] = sort(rej,'descend');
fprintf('  most-rejected pairs:');
for q = ord(1:min(5,numel(ord)))
  if rej(q) == 0, break; end
  fprintf(' %d-%d(%d)', P(q,1), P(q,2), rej(q));
end
fprintf('\n');

%% Fit the per-pass displacement against time and tide
tday = (ptime - min(ptime))/86400;
tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);

%% The same quantity from the single-reference pairing, for comparison
[Pm, Dm, ~, ~, ~] = load_network(ref_dir, PASS_NAME, REF_DEPTH, MAX_BASELINE);
Am = [];
if ~isempty(Pm)
  refp = mode(Pm(:,1));
  sec  = Pm(:,2).';
  Sm = nan(size(Dm,1), numel(sec));
  Sm(:, :) = Dm;
  Am = vdef.fitTideAdmittance(Sm, tday(sec)-min(tday(sec)), tide(sec));
end

span = max(tday)-min(tday);
% the single-reference comparison may be absent; [] makes prn print '-'
adm_m = []; adm_s = []; sec_m = []; sec_s = [];
if ~isempty(Am)
  adm_m = 1e3*REF_DEPTH*Am.admittance;
  adm_s = 1e3*REF_DEPTH*Am.admittance_std;
  sec_m = 1e3*REF_DEPTH*Am.trend*span;
  sec_s = 1e3*REF_DEPTH*Am.trend_std*span;
end
fprintf('\n%-28s %14s %14s\n','quantity (top 100 m)','network','single ref');
fprintf('%-28s %14d %14d\n','pairs used', round(median(N.n_used)), size(Pm,1));
prn('tide response [mm/m]', 1e3*REF_DEPTH*A.admittance, adm_m);
prn('  its 1-sigma [mm/m]', 1e3*REF_DEPTH*A.admittance_std, adm_s);
prn('secular over window [mm]', 1e3*REF_DEPTH*A.trend*span, sec_m);
prn('  its 1-sigma [mm]', 1e3*REF_DEPTH*A.trend_std*span, sec_s);

%% Is the limiting error per-PAIR or per-PASS?
% This is the question the network can answer and a single reference
% cannot. Redundancy averages down error that is independent per pair. It
% does NOTHING against error attached to a PASS: if epoch k is off by e_k
% then every pair touching it is off by +/-e_k, the network solves
% x_k = truth + e_k exactly, and no amount of extra pairs helps.
%
% The signature is whether the post-fit residual pattern over epochs
% REPEATS from block to block. Independent pair noise gives an
% uncorrelated pattern; a per-pass error gives the same passes off in the
% same direction everywhere along the line.
res_ep = nan(size(N.x));
for b = 1:size(N.x,1)
  y = N.x(b,:); ok = isfinite(y) & isfinite(tday) & isfinite(tide);
  if nnz(ok) < 6, continue; end
  X = [ones(nnz(ok),1), tday(ok).'-mean(tday(ok)), tide(ok).'-mean(tide(ok))];
  res_ep(b,ok) = (y(ok).' - X*(X\y(ok).')).';
end
good = find(sum(isfinite(res_ep),2) >= 6);
rho = [];
for a1 = 1:numel(good)
  for a2 = a1+1:numel(good)
    u = res_ep(good(a1),:); v = res_ep(good(a2),:);
    ok = isfinite(u) & isfinite(v);
    if nnz(ok) >= 6 && std(u(ok))>0 && std(v(ok))>0
      m = corrcoef(u(ok), v(ok)); rho(end+1) = m(1,2); %#ok<AGROW>
    end
  end
end
fprintf('\n===== is the error per-pair or per-pass? =====\n');
fprintf('block-to-block correlation of the per-epoch residual pattern:\n');
fprintf('  median %.2f over %d block pairs (0 = independent pair noise,\n', ...
  median(rho), numel(rho));
fprintf('  approaching 1 = the same passes are off everywhere, i.e. per-PASS)\n');
rms_ep = sqrt(mean(res_ep.^2,1,'omitnan'));
[~, wo] = sort(rms_ep,'descend');
fprintf('worst passes by residual [mm]:');
for q = wo(1:min(4,numel(wo)))
  if ~isfinite(rms_ep(q)), continue; end
  fprintf(' pass %d (%.2f)', q, 1e3*REF_DEPTH*rms_ep(q));
end
fprintf('\n');

%% Is the model too simple? Allow the tide response a PHASE LAG.
% Per-pair noise is only ~0.4 mm (the closure number) but the a+b*t+c*tide
% fit scatters at ~1.5 mm, and that extra is not common across blocks. One
% candidate is the model rather than the data: ice flexure is viscoelastic,
% so the strain response need not be IN PHASE with the tide, and a single
% c*tide term cannot represent a lag. Adding the quadrature component -
% d(tide)/dt, exactly 90 degrees out of phase for a sinusoid - lets the
% response take any phase.
tide_csv = fullfile(fileparts(mfilename('fullpath')),'diagnostics', ...
  'cats2008_apres_window.csv');
if exist(tide_csv,'file')
  % columns are datenum,iso,tide_m - the middle one is text, which
  % importdata cannot cope with, so read it explicitly
  fid = fopen(tide_csv,'r'); fgetl(fid);
  C = textscan(fid, '%f%s%f', 'Delimiter', ',');
  fclose(fid);
  tdn = C{1}; thh = C{3};
  assert(numel(tdn) == numel(thh) && numel(tdn) > 10, ...
    'failed to parse %s (%d/%d rows)', tide_csv, numel(tdn), numel(thh));
  % pass times as datenums, and the tide rate by central difference
  pdn = 719529 + ptime/86400;
  h_p  = interp1(tdn, thh, pdn, 'linear', NaN);
  dh_p = interp1(tdn(2:end-1), (thh(3:end)-thh(1:end-2))./(tdn(3:end)-tdn(1:end-2)), ...
    pdn, 'linear', NaN);
  r_in = nan(size(N.x,1),1); r_lag = nan(size(N.x,1),1);
  for b = 1:size(N.x,1)
    y = N.x(b,:);
    ok = isfinite(y) & isfinite(tday) & isfinite(h_p) & isfinite(dh_p);
    if nnz(ok) < 7, continue; end
    t0 = tday(ok).'-mean(tday(ok)); h0 = h_p(ok).'-mean(h_p(ok));
    q0 = dh_p(ok).'-mean(dh_p(ok)); yv = y(ok).';
    X1 = [ones(nnz(ok),1), t0, h0];
    X2 = [X1, q0];
    r_in(b)  = std(yv - X1*(X1\yv));
    r_lag(b) = std(yv - X2*(X2\yv));
  end
  ok = isfinite(r_in) & isfinite(r_lag);
  fprintf('\n===== does a tidal PHASE LAG explain the extra scatter? =====\n');
  fprintf('residual scatter, in-phase only : %.2f mm\n', 1e3*REF_DEPTH*median(r_in(ok)));
  fprintf('residual scatter, lag allowed   : %.2f mm\n', 1e3*REF_DEPTH*median(r_lag(ok)));
  fprintf('reduction %.0f%% over %d blocks\n', ...
    100*(1-median(r_lag(ok))/median(r_in(ok))), nnz(ok));
  fprintf(['A large reduction would mean the response is out of phase with\n' ...
    'the tide and the in-phase model was mis-specified; a small one means\n' ...
    'the scatter is data, not model.\n']);
else
  fprintf('\n(no CATS2008 csv at %s - skipping the phase-lag test)\n', tide_csv);
end

fprintf(['\nReading: the network carries a closure number a single reference\n' ...
  'cannot produce, and it names its own bad pairs. If it does NOT tighten\n' ...
  'the 1-sigma, that is the finding - more pairs cannot fix an error that\n' ...
  'belongs to the passes rather than to the pairings between them.\n']);

%% ========================================================================
function prn(lbl, a, b)
a = a(isfinite(a));
if isempty(b)
  fprintf('%-28s %14.2f %14s\n', lbl, median(a), '-');
else
  b = b(isfinite(b));
  fprintf('%-28s %14.2f %14.2f\n', lbl, median(a), median(b));
end
end

%% ========================================================================
function [P, D, W, along, nskip] = load_network(dirn, pn, REF_DEPTH, MAX_BASELINE)
P = []; D = []; W = []; along = []; nskip = 0;
f = dir(fullfile(dirn, [pn '_vvel_*.mat']));
for q = 1:numel(f)
  tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
    'tokens','once');
  if isempty(tok), continue; end
  o = load(fullfile(dirn, f(q).name));
  % a pair whose coalignment was rejected is NOT usable: it still carries
  % the multipass misalignment, which leaks into dtau at about 6%
  if ~vdef.pairAligned(o), nskip = nskip+1; continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, nskip = nskip+1; continue; end
  Nblk = numel(o.S1); sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH, 'linear', NaN)/REF_DEPTH;
  end
  if all(~isfinite(sv)), nskip = nskip+1; continue; end
  if isempty(D), D = sv; along = o.Along_track(:); else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'), 1e-3); %#ok<AGROW>
end
end
