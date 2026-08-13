%CLOSURE Does the strain measurement close around a loop?
%
%   THIS IS THE TRUTH-FREE QC FOR THE SECULAR TERM, the counterpart to
%   artefact_vs_signal.m for the tidal one. Run it before believing any
%   secular rate.
%
%   Results on 2026-08-13, with coalign_min_quality still at 0.50:
%     GL3  closure 0.96 mm rms, worst 1.69 mm
%     GL4  closure 3.23 mm rms, worst 7.05 mm
%   against an expected secular signal of about 2.4 mm over the window -
%   i.e. the internal inconsistency was as large as the signal. Raising
%   the quality floor to 0.85 cut GL4's worst case to 4.37 mm. Closure
%   quality also tracks the calibration split: GL3 carries a nonzero
%   coregistration_time_shift and closes 3x better than GL4, which does
%   not.
%
% THE TEST. Deformation is additive: going from pass i to pass k directly
% must equal going i->i+1->...->k step by step. The main pairing measures
% the direct long-baseline path (up to 3 days); the sequential pairing
% measures every short step (hours). Their difference is a CLOSURE
% RESIDUAL, and it needs no truth, no model and no tide - it is a pure
% internal-consistency check on the measurement.
%
% Why it matters here: closure residual is the standard signature of a
% phase-unwrapping or decorrelation bias in InSAR. Short pairs are highly
% coherent and unwrap cleanly; long pairs decorrelate. If the long pairs
% are biased, closure fails, and because the main pairing is ALL long
% pairs referenced to one epoch, that bias lands squarely in the secular
% term - which is exactly the term the five products disagree on.
%
% Only GL3 and GL4 have both pairings, which is enough: they are the two
% cleanest products (small cross-track baselines throughout).

vvel_dir  = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel_v3';
seq_dir   = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_vvel_seq_v3';
mp_dir    = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
REF_DEPTH = 100;
names = {'EAGER_2022_GL3','EAGER_2022_GL4'};

for n = 1:numel(names)
  pn = names{n};
  fprintf('\n===== %s =====\n', pn);

  [seq_i, seq_j, seq_s, seq_dt, seq_coh] = load_pairs(seq_dir, pn, REF_DEPTH);
  [mn_i,  mn_j,  mn_s,  mn_dt,  mn_coh ] = load_pairs(vvel_dir, pn, REF_DEPTH);
  if isempty(seq_i) || isempty(mn_i), fprintf('  missing a pairing\n'); continue; end
  ref = mode(mn_i);

  % cumulative strain from the sequential chain, relative to the reference
  Np = max([seq_i seq_j mn_i mn_j]);
  step = nan(1,Np-1);          % step(k) = strain from pass k to pass k+1
  for q = 1:numel(seq_i)
    if seq_j(q) == seq_i(q)+1, step(seq_i(q)) = seq_s(q); end
  end
  cum = nan(1,Np);             % cum(k) = strain from ref to pass k
  cum(ref) = 0;
  for k = ref+1:Np
    if k-1 >= 1 && isfinite(cum(k-1)) && isfinite(step(k-1)), cum(k) = cum(k-1) + step(k-1); end
  end
  for k = ref-1:-1:1
    if isfinite(cum(k+1)) && isfinite(step(k)), cum(k) = cum(k+1) - step(k); end
  end

  fprintf('  reference pass %d, %d sequential steps, %d direct pairs\n', ...
    ref, nnz(isfinite(step)), numel(mn_i));
  fprintf('%6s %9s %12s %12s %12s %9s %9s\n', 'pass','dt[d]', ...
    'direct[mm]','chained[mm]','closure[mm]','coh_dir','n_steps');
  cr = []; dtv = [];
  for q = 1:numel(mn_i)
    k = mn_j(q);
    if ~isfinite(cum(k)), continue; end
    d_mm = mn_s(q)*REF_DEPTH*1e3;
    c_mm = cum(k)*REF_DEPTH*1e3;
    nst = abs(k-ref);
    fprintf('%6d %9.2f %12.3f %12.3f %12.3f %9.3f %9d\n', ...
      k, mn_dt(q), d_mm, c_mm, d_mm-c_mm, mn_coh(q), nst);
    cr(end+1) = d_mm - c_mm; dtv(end+1) = mn_dt(q); %#ok<AGROW>
  end
  if numel(cr) >= 4
    ok = isfinite(cr) & isfinite(dtv);
    rr = corrcoef(abs(dtv(ok)), cr(ok));
    p = polyfit(dtv(ok), cr(ok), 1);
    fprintf('\n  closure residual: mean %+.2f mm, rms %.2f mm, max |%.2f| mm\n', ...
      mean(cr(ok)), sqrt(mean(cr(ok).^2)), max(abs(cr(ok))));
    fprintf('  vs |dt|: r = %+.2f ;  slope vs signed dt = %+.3f mm/day\n', ...
      rr(1,2), p(1));
    fprintf('  (a nonzero slope means the direct long-baseline pairs carry a\n');
    fprintf('   drift the short chained ones do not - i.e. a secular bias)\n');
  end
end

%% ========================================================================
function [pi_, pj_, s_, dt_, coh_] = load_pairs(dirn, pn, REF_DEPTH)
pi_=[]; pj_=[]; s_=[]; dt_=[]; coh_=[];
f = dir(fullfile(dirn, [pn '_vvel_*.mat']));
for i = 1:numel(f)
  tok = regexp(f(i).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], 'tokens','once');
  if isempty(tok), continue; end
  o = load(fullfile(dirn, f(i).name), 'depth_blk','dh_blk','S1','delta_t_sec','coh_blk');
  Nblk = numel(o.S1); sv = nan(1,Nblk);
  for b = 1:Nblk
    d = o.depth_blk(:,b); ok = isfinite(d) & isfinite(o.dh_blk(:,b));
    if ~any(ok) || max(d(ok)) < REF_DEPTH, continue; end
    sv(b) = interp1(d(ok), o.dh_blk(ok,b), REF_DEPTH, 'linear', NaN)/REF_DEPTH;
  end
  pi_(end+1) = str2double(tok{1}); pj_(end+1) = str2double(tok{2}); %#ok<AGROW>
  s_(end+1)  = median(sv(isfinite(sv))); %#ok<AGROW>
  dt_(end+1) = o.delta_t_sec/86400; %#ok<AGROW>
  coh_(end+1)= mean(o.coh_blk(:),'omitnan'); %#ok<AGROW>
end
end
