%FAR_END_CHECK How much weight will the last along-track bin carry?
%
%   WHY THIS EXISTS. Locating the ApRES sites (19 Aug 2026) put GA04 at
%   4.70 km along track, in the OUTERMOST stack bin, where the radar
%   reports +2.22 +/- 0.86 mm per metre of tide against ApRES -1.24 +/-
%   0.04. That is a ~4 sigma, opposite-sign disagreement, and it is about
%   to be drawn in the convergence figure - so the bin carrying it has to
%   be shown to be worth believing first. The far end is exactly where
%   this dataset is weakest: coherence thins with depth there, fewer lines
%   reach that far, and the GPS-curvature prediction is deliberately
%   trimmed 500 m from each end because a polynomial second derivative is
%   least constrained at the edges.
%
%   WHAT IT TESTS, all on the network products (CSARP_vvel_net):
%     1. WHICH LINES reach past 4 km, and what each says there - a bin
%        carried by one line is a single measurement, not a stack.
%     2. BLOCK QUALITY against the mid-line blocks: mean coherence, valid
%        depth samples, fitted depth span. If the far blocks are visibly
%        thinner, the formal sigma understates the real error.
%     3. BIN WIDTH: the same stack at 500 m and 1 km. A value that moves
%        a lot when neighbouring blocks are folded in is not resolved.
%        4.70 km is only bracketed by bin centres at 500 m (4.25 and
%        4.75); at 1 km the last centre is 4.50, so the site value there
%        is the 4.0-5.0 km bin that CONTAINS it, not an interpolation.
%        Every reported value says which of the two it is.
%     4. JACKKNIFE: drop each line in turn and restack. If one line
%        carries the sign, the +2.2 is that line's result, not the array's.
%
%   Run on the server:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('.../far_end_check.m')"

addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % +vdef

root    = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
mp_dir  = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
net_dir = fullfile(root,'CSARP_vvel_net');
PASS_NAMES = {'EAGER_2022_GL1','EAGER_2022_GL2', ...
              'EAGER_2022_GL3','EAGER_2022_GL4'};
REF_DEPTH = 100; MAX_BASELINE = 10;
FAR_KM = 4.0;                  % "far end" starts here
APRES_MM = -1.24; APRES_SE = 0.04; APRES_ALONG_KM = 4.70;

R = [];
for n = 1:numel(PASS_NAMES)
  S = one_line(PASS_NAMES{n}, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE);
  if isempty(S), continue; end
  if isempty(R), R = S; else, R(end+1) = S; end %#ok<AGROW>
end
assert(~isempty(R), 'no lines loaded');

%% 1. Which lines reach the far end, and what do they say there
fprintf('\n=== 1. per-line blocks beyond %.1f km ===\n', FAR_KM);
fprintf('%-6s %8s %9s %9s %8s %8s %8s %6s\n', ...
  'line','along_km','adm','adm_std','meancoh','nvalid','span_m','npair');
n_far_lines = 0;
for i = 1:numel(R)
  sel = find(R(i).along/1e3 >= FAR_KM & isfinite(R(i).adm));
  if ~isempty(sel), n_far_lines = n_far_lines + 1; end
  for b = sel(:).'
    fprintf('%-6s %8.2f %9.2f %9.2f %8.2f %8.0f %8.0f %6.0f\n', ...
      strrep(R(i).name,'EAGER_2022_',''), R(i).along(b)/1e3, ...
      R(i).adm(b), R(i).adm_std(b), R(i).coh(b), R(i).nvalid(b), ...
      R(i).span(b), R(i).npair(b));
  end
end
fprintf('%d of %d lines reach beyond %.1f km\n', n_far_lines, numel(R), FAR_KM);

%% 2. Far-end block quality against the mid-line blocks
fprintf('\n=== 2. block quality, far end vs mid line ===\n');
coh_f=[]; nv_f=[]; sp_f=[]; sg_f=[]; coh_m=[]; nv_m=[]; sp_m=[]; sg_m=[];
for i = 1:numel(R)
  a = R(i).along/1e3; ok = isfinite(R(i).adm);
  f = ok & a >= FAR_KM; m = ok & a > 1.0 & a < 3.5;
  coh_f=[coh_f;R(i).coh(f)]; nv_f=[nv_f;R(i).nvalid(f)]; sp_f=[sp_f;R(i).span(f)]; sg_f=[sg_f;R(i).adm_std(f)]; %#ok<AGROW>
  coh_m=[coh_m;R(i).coh(m)]; nv_m=[nv_m;R(i).nvalid(m)]; sp_m=[sp_m;R(i).span(m)]; sg_m=[sg_m;R(i).adm_std(m)]; %#ok<AGROW>
end
fprintf('%-10s %8s %8s\n','metric','far','mid');
fprintf('%-10s %8.2f %8.2f\n','mean coh', mean(coh_f,'omitnan'), mean(coh_m,'omitnan'));
fprintf('%-10s %8.0f %8.0f\n','n valid',  mean(nv_f,'omitnan'),  mean(nv_m,'omitnan'));
fprintf('%-10s %8.0f %8.0f\n','span m',   mean(sp_f,'omitnan'),  mean(sp_m,'omitnan'));
fprintf('%-10s %8.2f %8.2f\n','adm sigma',mean(sg_f,'omitnan'),  mean(sg_m,'omitnan'));
fprintf('n blocks: far %d, mid %d\n', nnz(isfinite(sg_f)), nnz(isfinite(sg_m)));

%% 3. Same stack at two bin widths
fprintf('\n=== 3. stack at 500 m vs 1 km bins, last 2 km ===\n');
for w = [0.5 1.0]
  [xc, sm, ss] = stack(R, w);
  keep = xc >= 3.0;
  fprintf('bin %.1f km:\n', w);
  kk = find(keep);
  for k = kk(:).'
    if isfinite(sm(k))
      fprintf('   %5.2f km  %+6.2f +/- %.2f\n', xc(k), sm(k), ss(k));
    end
  end
  [v, e, how] = at_site(xc, sm, ss, APRES_ALONG_KM, w);
  fprintf('  at the ApRES site (%.2f km): %+.2f +/- %.2f  -> %.1f sigma from ApRES %.2f  [%s]\n', ...
    APRES_ALONG_KM, v, e, abs(v-APRES_MM)/hypot(e,APRES_SE), APRES_MM, how);
end

%% 4. Jackknife: does one line carry the far-end sign?
fprintf('\n=== 4. jackknife of the value at the ApRES site ===\n');
[xc, sm, ss] = stack(R, 0.5);
fprintf('%-14s %+8s %8s  %s\n','dropped','value','sigma','from');
[v, e, how] = at_site(xc, sm, ss, APRES_ALONG_KM, 0.5);
fprintf('%-14s %+8.2f %8.2f  %s\n','none (all 4)', v, e, how);
for i = 1:numel(R)
  [xj, sj, sjs] = stack(R([1:i-1 i+1:end]), 0.5);
  [v, e, how] = at_site(xj, sj, sjs, APRES_ALONG_KM, 0.5);
  fprintf('%-14s %+8.2f %8.2f  %s\n', strrep(R(i).name,'EAGER_2022_',''), ...
    v, e, how);
end

%% ========================================================================
function [v, e, how] = at_site(xc, sm, ss, x0, w)
% Stacked value at along-track position x0, and a string saying where it
% came from. x0 = 4.70 km is only BRACKETED by bin centres at some widths:
% at w = 0.5 the centres 4.25 and 4.75 straddle it, but at w = 1.0 the last
% centre is 4.50, and interpolating there would be extrapolation. Rather
% than emit a silent NaN, fall back to the bin that CONTAINS x0 - it is a
% coarser answer to the same question, and the caller reports which it got.
ok = isfinite(sm) & isfinite(ss);
v = NaN; e = NaN; how = 'no stacked bin at this position';
if ~any(ok), return; end
xo = xc(ok); so = sm(ok); eo = ss(ok);
if numel(xo) >= 2 && x0 >= min(xo) && x0 <= max(xo)
  v = interp1(xo, so, x0, 'linear');
  e = interp1(xo, eo, x0, 'linear');
  how = 'interpolated between bin centres';
  return;
end
[dmin, k] = min(abs(xo - x0));
if dmin <= w/2
  v = so(k); e = eo(k);
  how = sprintf('%.2f-%.2f km bin', xo(k)-w/2, xo(k)+w/2);
end
end

%% ========================================================================
function [xc, sm, ss] = stack(R, w)
% Inverse-variance stack of the per-block admittance in bins of w km.
xmax = 0;
for i = 1:numel(R), xmax = max(xmax, max(R(i).along)/1e3); end
edges = 0:w:ceil(xmax/w)*w;
xc = edges(1:end-1) + w/2;
sm = nan(size(xc)); ss = nan(size(xc));
for k = 1:numel(xc)
  v = []; q = [];
  for i = 1:numel(R)
    al = R(i).along(:)/1e3; av = R(i).adm(:); as_ = R(i).adm_std(:);
    inb = al >= edges(k) & al < edges(k+1) & isfinite(av) & isfinite(as_) & as_ > 0;
    v = [v av(inb).']; q = [q 1./as_(inb).'.^2]; %#ok<AGROW>
  end
  if numel(v) >= 2
    sm(k) = sum(q.*v)/sum(q); ss(k) = sqrt(1/sum(q));
  end
end
end

%% ========================================================================
function S = one_line(pn, net_dir, mp_dir, REF_DEPTH, MAX_BASELINE)
% Per-block tidal admittance from the network product, plus the per-block
% quality metrics this diagnostic needs (coherence, valid-sample count and
% fitted depth span) that the figure script does not carry.
S = [];
f = dir(fullfile(net_dir, [pn '_vvel_*.mat']));
keep = ~cellfun('isempty', regexp({f.name}, ...
  ['^' regexptranslate('escape',pn) '_vvel_\d+_\d+\.mat$'], 'once'));
f = f(keep);
if isempty(f), return; end
L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
Np = numel(L.pass); elev = nan(1,Np); ptime = nan(1,Np);
for k = 1:Np
  elev(k)  = mean(L.pass(k).elev,'omitnan');
  ptime(k) = mean(L.pass(k).gps_time,'omitnan');
end

P = []; D = []; W = []; along = []; Nblk = 0;
coh = []; nvalid = []; span = [];
for q = 1:numel(f)
  tok = regexp(f(q).name, '_vvel_(\d+)_(\d+)\.mat$', 'tokens','once');
  o = load(fullfile(net_dir, f(q).name));
  if isfield(o,'coalign_applied') && ~o.coalign_applied, continue; end
  if max(abs(o.baseline_y)) > MAX_BASELINE, continue; end
  if Nblk == 0
    Nblk = numel(o.S1); along = o.Along_track(:);
    coh = zeros(Nblk,1); nvalid = zeros(Nblk,1); span = zeros(Nblk,1);
    npair = zeros(Nblk,1);
  end
  sv = nan(Nblk,1);
  for b = 1:Nblk
    d = o.depth_blk(:,b); okd = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(okd) || max(d(okd)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(okd), o.dh_blk(okd,b), REF_DEPTH,'linear',NaN)/REF_DEPTH;
  end
  if all(~isfinite(sv)), continue; end
  % accumulate quality over EVERY pair used, not just the first file -
  % a single pair is not representative of the block. Below the sv check,
  % so these averages describe the same pair set the admittance is
  % computed from
  for b = 1:Nblk
    d = o.depth_blk(:,b); okd = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(okd), continue; end
    coh(b) = coh(b) + mean(o.coh_blk(okd,b),'omitnan');
    nvalid(b) = nvalid(b) + nnz(okd);
    span(b) = span(b) + (max(d(okd)) - min(d(okd)));
    npair(b) = npair(b) + 1;
  end
  if isempty(D), D = sv; else, D(:,end+1) = sv; end %#ok<AGROW>
  P(end+1,:) = [str2double(tok{1}), str2double(tok{2})]; %#ok<AGROW>
  W(end+1) = max(mean(o.coh_blk(:),'omitnan'),1e-3); %#ok<AGROW>
end
if size(P,1) < 10, return; end

npair(npair == 0) = NaN;          % blocks no pair could measure
coh = coh ./ npair; nvalid = nvalid ./ npair; span = span ./ npair;

N = vdef.invertNetwork(P, D, struct('n_sigma',3,'weights',W,'n_epoch',Np));
tday = (ptime - min(ptime))/86400; tide = elev - mean(elev);
A = vdef.fitTideAdmittance(N.x, tday, tide);
S = struct('name',pn,'along',along, ...
  'adm', 1e3*REF_DEPTH*A.admittance(:), ...
  'adm_std', 1e3*REF_DEPTH*A.admittance_std(:), ...
  'coh', coh, 'nvalid', nvalid, 'span', span, 'npair', npair);
end
