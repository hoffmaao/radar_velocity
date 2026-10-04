function P = tidal_phase(opts)
%TIDAL_PHASE Does the column respond in phase with the tide? The lag from the quadrature fit.
%
%   An elastic column strains in phase with the tide, so its strain rate
%   leads the tide by exactly a quarter period. scripts/diagnostics/
%   tidal_stack.m fits every block and depth as Phi = kz*(a*dtide + b*dq +
%   v*dt), dq the quadrature tide (the tide rate times T/(2*pi)); a
%   response lagging the tide by tau has a + i*b = A*exp(-i*w*tau), w =
%   2*pi/T. This reads the lag off those cells.
%
%   ONE LAG FOR MANY CELLS. Each cell's (a, b) lies on the line through
%   the origin at angle -w*tau, whatever the sign and size of its A (the
%   response changes sign across the hinge, so angles of single cells
%   cannot simply be averaged). The common lag minimises the weighted
%   squared distance of the cells from that line,
%
%       chi2(tau) = sum (a sin(w tau) + b cos(w tau))^2 / s^2,
%       s^2 = (sd_a^2 + sd_b^2)/2,
%
%   and its 1-sigma interval is where chi2 rises by 1 (scaled by the
%   reduced chi2 at the minimum when that exceeds 1). Cells with no
%   response carry no information on the angle and only add a constant.
%   The quadrature amplitude alone is also tested against noise:
%   chi2_b = sum (b/sd_b)^2 over the cells, against its count.
%
%   WHICH INTERVAL TO QUOTE. The cells are not independent: every depth of
%   a block, and neighbouring blocks, are built from the same passes, so a
%   pass's error moves many cells together and the chi2 interval above is
%   too narrow (on the 3 Oct run per-line lags disagreed by several of
%   their own sigmas). The four lines share no passes, so they are four
%   independent replicates: the lag to quote is their mean and its standard
%   error from the line-to-line scatter (P.across).
%
%   THE CAVEAT. The tide at the pass times is a function of the hour of day
%   (R^2 = 0.99), so a lag is indistinguishable from a diurnal process that
%   peaks at another hour. A zero lag within error is the clean outcome.
%
%   opts: .files (tidal_stack_nozc.mat, tidal_stack_nozc_ref60.mat),
%   .labels, .zmin per file ([30 70], below each reference), .period_h
%   (24.84, as pass_tide), .near_km [-0.6 1.5] (where the bending is),
%   .out_dir. Returns P(file).line(k) and P(file).pooled with fields tau_h,
%   tau_lo, tau_hi, chi2r, chi2b, n, and the cells used, and P(file).across
%   (.tau_h, .se_h: the across-line mean lag and its standard error).

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath')); addpath(fileparts(fileparts(here)));
def = struct('out_dir', vdef.figureDir(), 'files', {{'tidal_stack_nozc.mat', 'tidal_stack_nozc_ref60.mat'}}, ...
  'labels', {{'reference 5 m', 'reference 60 m'}}, 'zmin', [30 70], 'period_h', 24.84, 'near_km', [-0.6 1.5]);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end
tg = linspace(-opts.period_h/4, opts.period_h/4, 1441);       % trial lags [h], +/- a quarter period
w = 2*pi/opts.period_h;
P = struct('label', {}, 'line', {}, 'pooled', {}, 'across', {}, 'near', {});
for f = 1:numel(opts.files)
  S = load(fullfile(opts.out_dir, opts.files{f})); OUT = S.OUT;
  assert(isfield(OUT, 'b_q'), '%s has no quadrature fit: rerun tidal_stack', opts.files{f});
  L = struct('name', {}, 'tau_h', {}, 'tau_lo', {}, 'tau_hi', {}, 'chi2r', {}, 'chi2b', {}, 'chi2a', {}, 'n', {}, 'r_qt', {});
  A = []; B = []; SA = []; SB = [];
  for k = 1:numel(OUT)
    O = OUT(k);
    sel = O.zsel(:) >= opts.zmin(f) & (O.x_sea(:).'/1e3 >= opts.near_km(1) & O.x_sea(:).'/1e3 <= opts.near_km(2));
    a = O.a_q(sel); b = O.b_q(sel); sa = O.a_q_sd(sel); sb = O.b_q_sd(sel);
    g = isfinite(a) & isfinite(b) & sa > 0 & sb > 0;
    a = a(g); b = b(g); sa = sa(g); sb = sb(g);
    e = fit_lag(a, b, sa, sb, tg, w);
    e.name = strrep(O.name, 'EAGER_2022_', ''); e.r_qt = O.r_qt;
    L(end+1) = orderfields(e, L); %#ok<AGROW>
    A = [A; a]; B = [B; b]; SA = [SA; sa]; SB = [SB; sb]; %#ok<AGROW>
  end
  pooled = fit_lag(A, B, SA, SB, tg, w);
  tl = [L.tau_h]; across = struct('tau_h', mean(tl), 'se_h', std(tl)/sqrt(numel(tl)));
  P(f) = struct('label', opts.labels{f}, 'line', L, 'pooled', pooled, 'across', across, ...
    'near', struct('a', A, 'b', B, 'sa', SA, 'sb', SB));
  fprintf('\n%s, near field (%.1f to %.1f km), below %d m; lag > 0 = response after the tide:\n', ...
    opts.labels{f}, opts.near_km, opts.zmin(f));
  fprintf('  %-6s %5s %7s %16s %7s %13s %13s\n', 'line', 'cells', 'r(q,T)', 'lag h (1 sigma)', 'chi2r', 'chi2_a/n', 'chi2_b/n');
  for k = 1:numel(L)
    fprintf('  %-6s %5d %+7.2f %+6.2f [%+5.2f %+5.2f] %7.2f %13.2f %13.2f\n', L(k).name, L(k).n, L(k).r_qt, ...
      L(k).tau_h, L(k).tau_lo, L(k).tau_hi, L(k).chi2r, L(k).chi2a, L(k).chi2b);
  end
  fprintf('  %-6s %5d %7s %+6.2f [%+5.2f %+5.2f] %7.2f %13.2f %13.2f\n', 'all', pooled.n, '', pooled.tau_h, ...
    pooled.tau_lo, pooled.tau_hi, pooled.chi2r, pooled.chi2a, pooled.chi2b);
  fprintf('  across the %d lines (independent passes): lag %+.2f +/- %.2f h (standard error)\n', ...
    numel(tl), across.tau_h, across.se_h);
end
fprintf(['\nchi2_a/n, chi2_b/n: in-phase and quadrature amplitudes against their jackknife sigmas (1 = noise).\n' ...
         'The per-line and pooled intervals treat cells as independent and are too narrow; quote the\n' ...
         'across-line value.\n' ...
         'A lag of 0 h means the strain rate leads the tide by exactly a quarter period (%.2f h).\n'], opts.period_h/4);
end

function e = fit_lag(a, b, sa, sb, tg, w)
s2 = (sa.^2 + sb.^2)/2;
c2 = arrayfun(@(t) sum((a*sin(w*t) + b*cos(w*t)).^2 ./ s2), tg);
[cmin, im] = min(c2); n = numel(a);
chi2r = cmin / max(n - 1, 1); sc = max(1, chi2r);
inside = c2 <= cmin + sc;
lo = find(~inside(1:im), 1, 'last'); hi = im - 1 + find(~inside(im:end), 1, 'first');
if isempty(lo), tlo = -Inf; else, tlo = tg(lo + 1); end    % reaches the edge: unbounded
if isempty(hi), thi = Inf; else, thi = tg(hi - 1); end
e = struct('tau_h', tg(im), 'tau_lo', tlo, 'tau_hi', thi, 'chi2r', chi2r, ...
  'chi2a', sum((a./sa).^2)/max(n, 1), 'chi2b', sum((b./sb).^2)/max(n, 1), 'n', n);
end
