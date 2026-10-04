function R = firn_modulus(opts)
%FIRN_MODULUS Where is the neutral plane, as a fraction of the local thickness?
%
%   In bending, plane sections stay plane, so horizontal strain is
%   kappa*(z - z_n) whatever the stiffness profile, and the neutral plane z_n
%   sits at the stiffness-weighted centroid of the column. Uniform stiffness
%   puts it at exactly half the LOCAL thickness, whatever that thickness is;
%   a firn softer than the ice pushes it deeper; a compliant lower column
%   pulls it shallower. Vertical strain follows -nu/(1-nu) times horizontal
%   strain, and below 60 m the density is within a few percent of ice, so
%   there the response relative to 60 m is
%
%       a(z) = A_b * [ (z^2 - 60^2)/2 - f*H_b*(z - 60) ],   70 m <= z < H_b - 15 m
%
%   with one amplitude per block, H_b the block's own ice thickness (the base
%   picked on its radargram, vdef.trackBase, stored by tidal_stack.m), and
%   the fraction f = z_n/H shared. The stack is masked below each block's
%   base, so the whole ice column enters wherever the ice is thick.
%
%   Uses the 60 m-referenced, trend-controlled stack with per-pass jackknife
%   sigmas, near-hinge blocks only; each line's far-field mean profile
%   (>= 1.5 km) is subtracted first, because the far field carries a
%   depth-shaped term of its own (opts.excess, default true).
%   Writes EAGER_2022_firn_modulus.png to vdef.figureDir.

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); addpath(here);
def = struct('out_dir', vdef.figureDir(), 'file', 'tidal_stack_nozc_ref60.mat', 'near_km', [-0.6 1.5], 'excess', true);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end, end
A = load(fullfile(opts.out_dir, opts.file)); OUT = A.OUT;
assert(isfield(OUT, 'base'), 'the stack predates the per-block ice base; rerun scripts/diagnostics/tidal_stack.m');
zr = 60; zmin = 70; f_grid = 0.33:0.002:0.75;

% ---- blocks: near-hinge excess over the line's far-field profile
blk = struct('z', {}, 'a', {}, 's', {}, 'H', {});
for k = 1:numel(OUT)
  O = OUT(k); z = O.zsel(:); zi = z >= zmin; xs = O.x_sea(:)/1e3; far = xs >= opts.near_km(2);
  af = 1e3*O.a_t(zi, far); sf = 1e3*O.a_t_sd(zi, far); wf = 1./sf.^2; wf(~isfinite(wf)) = 0;
  fm = sum(wf.*af, 2, 'omitnan')./sum(wf, 2); fs = 1./sqrt(sum(wf, 2));
  for b = find(xs >= opts.near_km(1) & xs < opts.near_km(2)).'
    a = 1e3*O.a_t(zi, b); s = 1e3*O.a_t_sd(zi, b);
    if opts.excess, a = a - fm; s = sqrt(s.^2 + fs.^2); end
    ok = isfinite(a) & isfinite(s) & s > 0; if nnz(ok) < 8, continue; end
    zz = z(zi); blk(end+1) = struct('z', zz(ok), 'a', a(ok), 's', s(ok), 'H', O.base(b)); %#ok<AGROW>
  end
end
shape = @(z, zn) (z.^2 - zr^2)/2 - zn*(z - zr);

% ---- profile chi2 over the fraction f
chi = zeros(size(f_grid)); nd = sum(arrayfun(@(q) numel(q.a), blk));
for g = 1:numel(f_grid)
  for q = 1:numel(blk)
    fz = shape(blk(q).z, f_grid(g)*blk(q).H); w = 1./blk(q).s.^2;
    Am = sum(w.*fz.*blk(q).a)/sum(w.*fz.^2); chi(g) = chi(g) + sum(((blk(q).a - Am*fz)./blk(q).s).^2);
  end
end
chi0 = sum(arrayfun(@(q) sum((q.a./q.s).^2), blk));
[cmin, im] = min(chi); f_best = f_grid(im);
red = cmin / max(nd - numel(blk) - 1, 1); dchi = (chi - cmin)/max(red, 1);
f1 = f_grid(dchi <= 1); f2 = f_grid(dchi <= 4);
Hm = median([blk.H]);
% soft firn: E ~ rho^3 on the Herron-Langway column, for the median thickness
par = vdef.defaultParams(); par.max_depth = Hm; par.nz = 2901; P = vdef.firnColumn(par);
f_firn = trapz(P.d, P.rho.^3 .* P.d)/trapz(P.d, P.rho.^3)/Hm;
R = struct('f_best', f_best, 'f_1sig', [min(f1) max(f1)], 'f_2sig', [min(f2) max(f2)], 'f_firn', f_firn, ...
  'H_median', Hm, 'H_range', [min([blk.H]) max([blk.H])], 'chi', cmin, 'red_chi', red, 'chi_nobend', chi0, ...
  'nblocks', numel(blk), 'lower_over_upper', (4*f_best - 1)/(3 - 4*f_best));
fprintf('near hinge: %d blocks, thickness %.0f-%.0f m (median %.0f)\n', R.nblocks, R.H_range, Hm);
fprintf('  z_n/H = %.3f (1 sigma %.3f-%.3f, 2 sigma %.3f-%.3f); uniform 0.500 (delta-chi2 %.1f); soft firn %.3f\n', ...
  f_best, R.f_1sig, R.f_2sig, interp1(f_grid, dchi, 0.5), f_firn);
fprintf('  = %.0f m at the median thickness; chi2 %.1f (reduced %.2f) with bending, %.1f with none\n', f_best*Hm, cmin, red, chi0);

% ---- figure
[h, GRL] = grl_figure(170, 92); set(0, 'CurrentFigure', h);
c_dat = [0.30 0.30 0.30]; c_fit = [0.698 0.094 0.169]; c_uni = [0.165 0.471 0.839]; c_firn = [0.85 0.55 0.10]; soft = [0.5 0.5 0.5];
% (a) every block's profile in depth-fraction, amplitude normalised, binned
ax1 = axes('parent', h, 'Position', [0.08 0.25 0.36 0.66]); hold(ax1, 'on');
zeta = []; v = []; sv = [];
for q = 1:numel(blk)
  zn = f_best*blk(q).H; fz = shape(blk(q).z, zn); w = 1./blk(q).s.^2; Am = sum(w.*fz.*blk(q).a)/sum(w.*fz.^2);
  S0 = abs(shape(blk(q).H - 15, zn));
  zeta = [zeta; blk(q).z/blk(q).H]; v = [v; blk(q).a/Am/S0]; sv = [sv; blk(q).s/abs(Am)/S0]; %#ok<AGROW>
end
edges = 0.22:0.04:0.96; zc = (edges(1:end-1) + edges(2:end))/2; ym = nan(size(zc)); ye = ym;
for i = 1:numel(zc)
  in = zeta >= edges(i) & zeta < edges(i+1) & isfinite(v); if ~any(in), continue; end
  w = 1./sv(in).^2; ym(i) = sum(w.*v(in))/sum(w); ye(i) = 1/sqrt(sum(w));
end
ok = isfinite(ym);
plot(ax1, [ym(ok) - ye(ok); ym(ok) + ye(ok)], [zc(ok); zc(ok)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hd = plot(ax1, ym(ok), zc(ok), 'o', 'Color', c_dat, 'MarkerFaceColor', c_dat, 'MarkerSize', 3.5);
zf = linspace(zr/Hm, 1 - 15/Hm, 200); zfm = zf*Hm; S0m = abs(shape(Hm - 15, f_best*Hm));
hb = plot(ax1, shape(zfm, f_best*Hm)/S0m, zf, '-', 'Color', c_fit, 'LineWidth', 1.6);
fu = shape(zc(ok)*Hm, 0.5*Hm)/S0m; ku = sum(fu.*ym(ok)./ye(ok).^2)/sum(fu.^2./ye(ok).^2);
hu = plot(ax1, ku*shape(zfm, 0.5*Hm)/S0m, zf, '--', 'Color', c_uni, 'LineWidth', 1.4);
plot(ax1, [-2 3], [f_best f_best], ':', 'Color', c_fit, 'LineWidth', 0.8);
plot(ax1, [-2 3], [0.5 0.5], ':', 'Color', c_uni, 'LineWidth', 0.8);
plot(ax1, [0 0], [0 1], '-', 'Color', soft, 'LineWidth', 0.5);
set(ax1, 'YDir', 'reverse', 'YLim', [0.2 1], 'XLim', [-2 3], 'FontSize', 7.5, 'Box', 'off', 'TickDir', 'out', 'Layer', 'top');
grid(ax1, 'on'); set(ax1, 'GridAlpha', 0.1);
xlabel(ax1, 'Response relative to 60 m (block amplitude = 1)', 'FontSize', 7.5);
ylabel(ax1, 'Depth / local ice thickness', 'FontSize', 7.5);
title(ax1, sprintf('(a)  Near-hinge excess, %d blocks, H = %.0f-%.0f m', R.nblocks, R.H_range), 'FontWeight', 'normal', 'FontSize', 8);
legend(ax1, [hd hb hu], {'data (binned)', sprintf('fit, z_n/H = %.2f', f_best), 'uniform, z_n/H = 0.50'}, ...
  'FontSize', 6.5, 'Box', 'off', 'Position', [0.29 0.62 0.10 0.13]);
% (b) fit score against the fraction
ax2 = axes('parent', h, 'Position', [0.58 0.25 0.38 0.66]); hold(ax2, 'on');
patch(ax2, [R.f_1sig fliplr(R.f_1sig)], [0 0 9 9], c_fit, 'FaceAlpha', 0.10, 'EdgeColor', 'none');
plot(ax2, f_grid, dchi, '-', 'Color', c_fit, 'LineWidth', 1.6);
plot(ax2, [0.5 0.5], [0 9], '--', 'Color', c_uni, 'LineWidth', 1.2);
text(ax2, 0.505, 8.5, {'uniform', 'stiffness'}, 'Color', c_uni, 'FontSize', 6.5, 'VerticalAlignment', 'top');
plot(ax2, f_firn*[1 1], [0 9], '--', 'Color', c_firn, 'LineWidth', 1.2);
text(ax2, f_firn + 0.005, 6.5, {'soft firn', 'E \propto \rho^3'}, 'Color', c_firn, 'FontSize', 6.5, 'VerticalAlignment', 'top');
text(ax2, 0.62, 2.6, {sprintf('measured z_n/H = %.2f', f_best), sprintf('(%.0f m of %.0f m)', f_best*Hm, Hm)}, ...
  'Color', c_fit, 'FontSize', 6.5, 'VerticalAlignment', 'bottom');
plot(ax2, f_grid([1 end]), [1 1], ':', 'Color', soft); plot(ax2, f_grid([1 end]), [4 4], ':', 'Color', soft);
set(ax2, 'XLim', f_grid([1 end]), 'YLim', [0 9], 'FontSize', 7.5, 'Box', 'off', 'TickDir', 'out');
grid(ax2, 'on'); set(ax2, 'GridAlpha', 0.1);
xlabel(ax2, 'Neutral plane, z_n / local ice thickness', 'FontSize', 7.5); ylabel(ax2, '\Delta\chi^2 (rescaled)', 'FontSize', 7.5);
title(ax2, '(b)  How well does each neutral plane fit?', 'FontWeight', 'normal', 'FontSize', 8);
annotation(h, 'textbox', [0.08 0.0 0.88 0.09], 'String', ...
  ['Trend-controlled coherent stack referenced at 60 m and masked below each block''s picked ice base. Near-hinge blocks (-0.6 to 1.5 km) minus ' ...
   'their line''s far-field profile, fitted with their own bending amplitude and one shared neutral plane at a fixed fraction of local thickness.'], ...
  'EdgeColor', 'none', 'FontSize', 6.5, 'Color', soft, 'VerticalAlignment', 'bottom');
out_fn = fullfile(opts.out_dir, 'EAGER_2022_firn_modulus.png');
print(h, out_fn, '-dpng', sprintf('-r%d', GRL.dpi)); close(h);
fprintf('Wrote %s\n', out_fn);
end
