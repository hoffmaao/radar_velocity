%TEST_BEAM_FLEXURE Unit test for vdef.beamFlexure and vdef.invertElasticModulus.
%
%   What the flexure chain has to make good: reproduce the textbook solution
%   for a beam it can be checked against analytically, converge at the order
%   the discretisation claims, recover a known modulus from a flexure
%   profile, and be honest about the two ways the inverse problem is not
%   determined - the E*-h^3 trade-off and a survey window that does not span
%   the whole flexure zone.
%
%   Covers:
%     1. the constant-thickness solution matches the analytic one, and the
%        four boundary conditions are satisfied
%     2. second-order convergence under grid refinement
%     3. a varying thickness produces a measurably different beam - the
%        whole reason for solving numerically rather than using (1)
%     4. round trip: a known E* is recovered from a well-sampled profile
%     5. the answer depends on the profile SHAPE alone, so an arbitrary
%        rescaling of the data leaves E* untouched. This is what lets a
%        tidal admittance normalised to a line mean be inverted at all.
%     6. the E*-h^3 trade-off is exactly that: assuming ice 10% thicker
%        moves E* by 1.1^-3, matching R.dlogE_dlogh
%     7. a partial window - flexure observed without either flat end, the
%        McMurdo case - still runs but reports a much wider interval, so
%        the ill-posedness shows up in the error bar rather than in a
%        confidently wrong number
%     8. observations that are block MEANS are modelled as block means.
%        Comparing them to the deflection at the block centre instead is a
%        systematic an order of magnitude above the formal error here
%     9. running x landward instead of seaward is refused, not fitted
%    10. joint fit with the englacial strain admittance: one shared
%        amplitude explains both observables of the same beam, and E* is
%        recovered with both per-dataset chi-squareds at order 1
%    11. the strain data carry AMPLITUDE information the normalised shape
%        cannot: on the partial window of check 7 the joint interval
%        closes to a fraction of the shape-only one. This is the reason to
%        include the radar at all, so it is asserted, not assumed
%    12. LOCAL E* (the E_patch mode) is honest about its resolution: a
%        patch in the bending zone recovers constant truth without
%        inventing structure, a genuinely stiffened patch is recovered,
%        and a patch out on the flat shelf - where moment and curvature
%        both vanish - comes back unconstrained rather than confidently
%        wrong. This is what makes mapping E* along the line defensible
%        where it is, and visibly meaningless where it is not
%
%   Runs in MATLAB or Octave:
%     docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work/scripts \
%       gnuoctave/octave:latest octave --no-gui test_beam_flexure.m

addpath(fileparts(fileparts(mfilename('fullpath'))));   % +vdef
rand('seed', 11); randn('seed', 11);   %#ok<RAND>

C     = vdef.constants();
NU    = 0.3;
K     = C.rho_sea * C.g;                 % foundation modulus [Pa/m]
rigid = @(E,h) E .* h.^3 / (12*(1 - NU^2));
flex  = @(E,h) (4*rigid(E,h) / K).^(1/4);           % flexural length [m]

%% 1. Analytic constant-thickness solution
% For uniform D on a semi-infinite beam clamped at x = 0, the solution of
% (D w'')'' + k w = k A0 is w = A0 [1 - exp(-x/l)(cos(x/l) + sin(x/l))],
% l = (4D/k)^(1/4). Taking the domain out to ~10 l makes the far-field
% boundary conditions the model imposes exact to a few parts in 10^4.
E0 = 1.0e9; h0 = 300; A0 = 0.85;
l0 = flex(E0, h0);
L  = 10*l0;
dx = 20;
x  = (0:dx:L).';
F  = vdef.beamFlexure(x, h0, E0, struct('A0', A0));

s   = x/l0;
wan = A0*(1 - exp(-s).*(cos(s) + sin(s)));
err = max(abs(F.w - wan))/A0;
fprintf('1. flexural length %.0f m, %d nodes, solve residual %.2e\n', ...
  l0, numel(x), F.residual);
fprintf('   max |numeric - analytic| = %.3e of A0\n', err);
assert(err < 2e-3, 'analytic mismatch %.3e is too large', err);
assert(F.residual < 1e-10, 'the linear system was not solved (residual %.2e)', F.residual);
assert(abs(F.w(1)) < 1e-12*A0, 'w(0) = 0 violated');
assert(abs(F.w(end) - A0) < 1e-12*A0, 'w(L) = A0 violated');
assert(abs(F.dwdx(1))*l0/A0 < 1e-10, 'dw/dx(0) = 0 violated');
assert(abs(F.dwdx(end))*l0/A0 < 1e-10, 'dw/dx(L) = 0 violated');

%% 2. Second-order convergence
e = zeros(1,3);
dxs = [80 40 20];
for q = 1:3
  xq = (0:dxs(q):L).';
  Fq = vdef.beamFlexure(xq, h0, E0, struct('A0', A0));
  sq = xq/l0;
  e(q) = max(abs(Fq.w - A0*(1 - exp(-sq).*(cos(sq) + sin(sq)))))/A0;
end
r = e(1:2)./e(2:3);
fprintf('2. errors %.2e %.2e %.2e, ratios %.2f %.2f (expect ~4)\n', e, r);
assert(all(r > 3.2) && all(r < 5.0), ...
  'convergence ratios %.2f %.2f are not second order', r);

%% 3. Variable thickness changes the beam
% Thickness enters as h^3, so a gradient the eye would call mild is not a
% small perturbation: this is why the paper solves a variable-D beam rather
% than picking a mean thickness.
hlin = 380 - 160*(x/L);                      % 380 m at the clamp to 220 m
Fv   = vdef.beamFlexure(x, hlin, E0, struct('A0', A0));
Fm   = vdef.beamFlexure(x, mean(hlin), E0, struct('A0', A0));
dmax = max(abs(Fv.w - Fm.w))/A0;
fprintf('3. h %.0f-%.0f m vs its mean: profiles differ by %.1f%% of A0\n', ...
  hlin(1), hlin(end), 100*dmax);
assert(dmax > 0.02, ...
  'a 380-220 m thickness gradient moved the beam by only %.3f of A0', dmax);

%% 4. Round trip on a well-sampled profile
% The ICESat-2 geometry: the track runs from flat grounded ice, through the
% flexure zone, out to flat floating ice, so both ends of the profile are
% observed and the clamp is inside the window.
%
% The truth comes from one forward solve on a 10 m grid, far finer than any
% the inversion builds for itself and running well past the last
% observation, so a round trip tests the fit rather than a grid coincidence.
E_true = 2.4e9; h_true = 280; x0_true = 3000;
xt = (x0_true : 10 : 26000).';
Ft = vdef.beamFlexure(xt, h_true, E_true, struct('A0', 1.0));

xo = (0:200:14000).';
W  = zeros(size(xo));
in = xo >= x0_true;
W(in) = interp1(Ft.x, Ft.w, xo(in), 'linear');
sig = 0.01;
wo = W + sig*randn(size(W));

iopts = struct('h', h_true, 'sigma', sig*ones(size(wo)), ...
               'E_grid', logspace(log10(0.2e9), log10(20e9), 41), ...
               'x0_grid', x0_true + (-2000:250:2000));
R = vdef.invertElasticModulus(xo, wo, iopts);
fprintf(['4. E* = %.2f GPa [%.2f, %.2f] (truth %.2f), x0 = %.0f m ' ...
         '(truth %.0f)\n'], R.E/1e9, R.E_lo/1e9, R.E_hi/1e9, E_true/1e9, ...
         R.x0, x0_true);
fprintf('   rms %.4f, chi2red %.2f, %.1f flexural lengths of pad, %d local minima\n', ...
  R.rms, R.chi2red, R.n_lambda, R.n_local_min);
assert(R.interior, 'the minimum ran to a grid edge');
assert(abs(log(R.E/E_true)) < log(1.15), ...
  'E* off by %.0f%%', 100*(R.E/E_true - 1));
% Two formal sigmas, not one. The reported interval IS one sigma, so a
% one-sigma coverage check on a single fixed noise realisation is a coin
% flip dressed up as an assertion - it would fail about a third of the time
% for a perfectly correct inversion.
nsig = abs(log(R.E/E_true)) / log(R.E_hi/R.E);
fprintf('   truth is %.1f formal sigma from the fit\n', nsig);
assert(nsig < 2, ...
  'the truth %.2f GPa is %.1f formal sigma from the fitted %.2f GPa', ...
  E_true/1e9, nsig, R.E/1e9);
assert(abs(R.x0 - x0_true) < 500, 'x0 off by %.0f m', R.x0 - x0_true);
assert(R.chi2red < 2.5, 'reduced chi-squared %.2f - the model does not fit', R.chi2red);
assert(R.n_local_min == 1, ...
  ['a clean synthetic reported %d rival minima; the rival filter is ' ...
   'letting the valley walls through and the flag is then worthless'], ...
  R.n_local_min);

%% 5. Only the shape is fitted
% The amplitude is eliminated in closed form, so multiplying the data by any
% constant must leave E* bit-identical. A GPS tidal admittance is normalised
% by a line mean rather than by the tide, and this is what makes it a
% legitimate input despite that.
Rs = vdef.invertElasticModulus(xo, 7.3*wo, ...
       struct('h', h_true, 'E_grid', iopts.E_grid, 'x0_grid', iopts.x0_grid));
Ru = vdef.invertElasticModulus(xo, wo, ...
       struct('h', h_true, 'E_grid', iopts.E_grid, 'x0_grid', iopts.x0_grid));
fprintf('5. data x7.3: E* %.4f -> %.4f GPa, amplitude %.3f -> %.3f\n', ...
  Ru.E/1e9, Rs.E/1e9, Ru.amplitude, Rs.amplitude);
assert(abs(Rs.E/Ru.E - 1) < 1e-12, 'rescaling the data moved E* by %.2e', ...
  Rs.E/Ru.E - 1);
assert(abs(Rs.amplitude/Ru.amplitude - 7.3) < 1e-9, ...
  'the fitted amplitude did not follow the data scale');

%% 6. The E*-h^3 trade-off, measured
% Flexure constrains D = E h^3/(12(1-nu^2)) and nothing finer, so assuming
% thicker ice must return a proportionally softer modulus. Noiseless data,
% so what is left is the trade-off itself and not fit scatter.
Rh1 = vdef.invertElasticModulus(xo, W, ...
        struct('h', h_true,     'E_grid', iopts.E_grid, 'x0_grid', iopts.x0_grid));
Rh2 = vdef.invertElasticModulus(xo, W, ...
        struct('h', 1.1*h_true, 'E_grid', iopts.E_grid, 'x0_grid', iopts.x0_grid));
ratio  = Rh2.E / Rh1.E;
expect = 1.1^Rh1.dlogE_dlogh;
fprintf('6. h x1.1: E* %.3f -> %.3f GPa, ratio %.4f (h^-3 predicts %.4f)\n', ...
  Rh1.E/1e9, Rh2.E/1e9, ratio, expect);
assert(abs(ratio/expect - 1) < 0.04, ...
  'ratio %.4f departs from the h^-3 prediction %.4f', ratio, expect);
assert(abs(Rh2.D/Rh1.D - 1) < 0.04, ...
  'the two fits should agree on rigidity, not modulus (D differs by %.1f%%)', ...
  100*(Rh2.D/Rh1.D - 1));

%% 7. A partial window is honest about it
% The McMurdo geometry: ~5 km of track inside the flexure zone, reaching
% neither the flat grounded end nor the flat floating one, so the clamp is
% outside the data and E* trades off against where it sits. The inversion
% should still run - and should say so in the width of its interval.
xp = (0:500:5000).';
Wp = interp1(Ft.x, Ft.w, xp + x0_true + 1500, 'linear');
Wp = Wp / mean(Wp);                          % a line-mean normalised a(x)
sig_p = 0.02;
wp = Wp + sig_p*randn(size(Wp));
Rp = vdef.invertElasticModulus(xp, wp, ...
       struct('h', h_true, 'sigma', sig_p*ones(size(wp)), ...
              'E_grid', logspace(log10(0.2e9), log10(20e9), 41), ...
              'x0_grid', -6000:250:0));
span_full = log10(R.E_hi/R.E_lo);
span_part = log10(Rp.E_hi/Rp.E_lo);
fprintf(['7. partial window: E* = %.2f GPa [%.2f, %.2f], x0 = %.0f m ' ...
         '(truth %.0f)\n'], Rp.E/1e9, Rp.E_lo/1e9, Rp.E_hi/1e9, Rp.x0, -1500);
fprintf('   interval spans %.2f decades against %.2f for the full profile\n', ...
  span_part, span_full);
% The bias is the point of this check. A window that sees neither flat end
% cannot separate E* from a clamp position it never observes, and the fit
% comes back tens of percent off. What it must not do is come back tens of
% percent off while claiming a few percent of precision - so the bias and
% the half-width are printed together, and the half-width has to cover it.
fprintf('   bias %+.0f%% of truth, against a formal half-width of %.0f%%\n', ...
  100*(Rp.E/E_true - 1), 100*(Rp.E_hi/Rp.E - 1));
assert(isfinite(Rp.E), 'the partial window returned no modulus at all');
nsig_p = abs(log(Rp.E/E_true)) / log(Rp.E_hi/Rp.E);
assert(nsig_p < 2, ...
  ['the partial window is %.1f formal sigma from the truth - it is being ' ...
   'confidently wrong rather than honestly uncertain'], nsig_p);
assert(isnan(span_part) || span_part > 1.5*span_full, ...
  ['the partial window reported an interval %.2f decades wide, no wider ' ...
   'than the full profile at %.2f - the ill-posedness is not being seen'], ...
  span_part, span_full);

%% 8. Block-averaged observations have to be modelled as block averages
% Every a(x) this project produces is a mean over a 500 m along-track block,
% and the deflection is curved on the scale of a flexural length, so the
% block mean is not the deflection at the block centre. Comparing it to the
% centre value is a systematic the fit cannot absorb - and it does not
% announce itself as a bad fit unless the formal errors are small enough to
% see it, which here they are.
BW = 500;
xblk = (3200 : BW : 12200).';
Wb = zeros(size(xblk));
% A brute-force reference average, fine enough that its own quadrature error
% is nine orders below the noise. The model's block average has to earn its
% agreement rather than share a rule with the thing checking it.
uu = linspace(-0.5, 0.5, 2001);
for q = 1:numel(xblk)
  xs = xblk(q) + BW*uu;
  ws = zeros(size(xs));
  inb = xs >= x0_true;
  ws(inb) = interp1(Ft.x, Ft.w, xs(inb), 'linear');
  Wb(q) = mean(ws);
end
sigb  = 1e-4;                                % small enough to see the bias
wb    = Wb + sigb*randn(size(Wb));
bopts = struct('h', h_true, 'sigma', sigb*ones(size(wb)), ...
               'E_grid', logspace(log10(0.5e9), log10(12e9), 41), ...
               'x0_grid', x0_true + (-1500:250:1500));
Rpt = vdef.invertElasticModulus(xblk, wb, bopts);
bopts.avg_width = BW;
Rbk = vdef.invertElasticModulus(xblk, wb, bopts);
fprintf(['8. %.0f m blocks: point-sampled E* = %.2f GPa (chi2red %.0f), ' ...
         'block-averaged %.2f GPa (chi2red %.1f), truth %.2f\n'], ...
  BW, Rpt.E/1e9, Rpt.chi2red, Rbk.E/1e9, Rbk.chi2red, E_true/1e9);
assert(Rbk.chi2red < 3, ...
  'block-averaged model still misfits at chi2red %.1f', Rbk.chi2red);
assert(Rpt.chi2red > 10*Rbk.chi2red, ...
  ['point sampling should misfit far worse than block averaging here ' ...
   '(%.1f vs %.1f); if it does not, this check has stopped testing ' ...
   'anything'], Rpt.chi2red, Rbk.chi2red);
assert(abs(log(Rbk.E/E_true)) < log(1.10), ...
  'block-averaged E* off by %.0f%%', 100*(Rbk.E/E_true - 1));

%% 9. A landward-running coordinate is refused
threw = false;
try
  vdef.invertElasticModulus(xo, flipud(wo), struct('h', h_true));
catch ME
  threw = true;
  assert(~isempty(strfind(upper(ME.message), 'SEAWARD')), ...
    'wrong error for a flipped profile: %s', ME.message);
end
assert(threw, 'a landward-running x was accepted instead of refused');
fprintf('9. a landward-running coordinate is refused\n');

%% 10. Joint fit with the englacial strain admittance
% The beam's internal bending strain, integrated over the top zr metres,
% is dh(zr) = A0 * K2 * w''(x) with K2 = nu/(1-nu)*(z_n*zr - zr^2/2). The
% check-4 shape data are in absolute metres (amplitude 1), so the strain
% observable here is K2*w'' directly; both truths come from the same fine
% beam Ft, block-averaged by the same brute-force rule as check 8.
ZR = 100; ZN = h_true/2;
K2 = (0.3/0.7) * (ZN*ZR - ZR^2/2);
xstr = (3250:500:13750).';
Wpp  = zeros(size(xstr));
for q = 1:numel(xstr)
  xsmp = xstr(q) + BW*uu;
  ws = zeros(size(xsmp));
  inb = xsmp >= x0_true;
  ws(inb) = interp1(Ft.x, Ft.d2wdx2, xsmp(inb), 'linear');
  Wpp(q) = mean(ws);
end
sig_s = 2e-4;                               % 0.2 mm, ApRES-class
dstr  = K2*Wpp + sig_s*randn(size(Wpp));
Rj = vdef.invertElasticModulus(xo, wo, ...
  struct('h', h_true, 'sigma', sig*ones(size(wo)), ...
         'E_grid', iopts.E_grid, 'x0_grid', iopts.x0_grid, ...
         'strain', struct('x', xstr, 'y', dstr, ...
                          'sigma', sig_s*ones(size(dstr)), ...
                          'ref_depth', ZR, 'avg_width', BW)));
fprintf(['10. joint: E* = %.2f GPa [%.2f, %.2f] (truth %.2f), ' ...
         'chi2 shape %.2f, strain %.2f (%d pts)\n'], ...
  Rj.E/1e9, Rj.E_lo/1e9, Rj.E_hi/1e9, E_true/1e9, ...
  Rj.chi2_shape, Rj.chi2_strain, Rj.n_strain);
assert(Rj.has_strain && Rj.n_strain == numel(xstr), 'strain dataset not carried');
assert(abs(log(Rj.E/E_true)) < log(1.12), 'joint E* off by %.0f%%', ...
  100*(Rj.E/E_true - 1));
nsig_j = abs(log(Rj.E/E_true)) / log(Rj.E_hi/Rj.E);
assert(nsig_j < 2, 'joint truth at %.1f sigma', nsig_j);
assert(Rj.chi2_shape < 2.5 && Rj.chi2_strain < 2.5, ...
  'a dataset misfits: shape %.2f, strain %.2f', Rj.chi2_shape, Rj.chi2_strain);

%% 11. Strain closes the partial window
% Check 7's geometry again - neither flat end observed - now with the
% strain admittance alongside the normalised shape. The shape was
% normalised by its own window mean, so the true shared amplitude is
% 1/mean(w) over the shape blocks, and the strain data are built with
% exactly that scale, the way a real line-mean-normalised a(x) and an
% absolute radar admittance relate.
amp_p = 1/mean(interp1(Ft.x, Ft.w, xp + x0_true + 1500, 'linear'));
wpp_p = interp1(Ft.x, Ft.d2wdx2, xp + x0_true + 1500, 'linear');
dstr_p = amp_p*K2*wpp_p + sig_s*randn(size(wpp_p));
Rpj = vdef.invertElasticModulus(xp, wp, ...
  struct('h', h_true, 'sigma', sig_p*ones(size(wp)), ...
         'E_grid', logspace(log10(0.2e9), log10(20e9), 41), ...
         'x0_grid', -6000:250:0, ...
         'strain', struct('x', xp, 'y', dstr_p, ...
                          'sigma', sig_s*ones(size(dstr_p)), 'ref_depth', ZR)));
span_joint = log10(Rpj.E_hi/Rpj.E_lo);
fprintf(['11. partial window joint: E* = %.2f GPa [%.2f, %.2f], x0 = %.0f m; ' ...
         'interval %.2f decades vs %.2f shape-only\n'], ...
  Rpj.E/1e9, Rpj.E_lo/1e9, Rpj.E_hi/1e9, Rpj.x0, span_joint, span_part);
nsig_pj = abs(log(Rpj.E/E_true)) / log(Rpj.E_hi/Rpj.E);
assert(nsig_pj < 2, 'joint partial-window truth at %.1f sigma', nsig_pj);
assert(span_joint < 0.6*span_part, ...
  ['the strain data did not tighten the partial window (%.2f vs %.2f ' ...
   'decades) - the amplitude constraint is not reaching the fit'], ...
  span_joint, span_part);

%% 12. Local E* is honest about where it is resolved
% The E_patch mode, on the full-window joint data of check 10. The clamp
% is FIXED at the global fit's answer for every patch - scanning patches
% with a free clamp lets x0 chase each patch along the trade-off valley.
popts = struct('h', h_true, 'sigma', sig*ones(size(wo)), ...
  'E_grid', logspace(log10(E_true/12), log10(E_true*12), 31), ...
  'x0_grid', Rj.x0, ...
  'strain', struct('x', xstr, 'y', dstr, 'sigma', sig_s*ones(size(dstr)), ...
                   'ref_depth', ZR, 'avg_width', BW));

% (a) constant truth, patch in the bending zone: no invented structure
popts.E_patch = struct('lo', 3000, 'hi', 5000, 'E_ref', Rj.E);
Pa = vdef.invertElasticModulus(xo, wo, popts);
fprintf('12a. bending-zone patch on constant truth: E* = %.2f [%.2f, %.2f] (truth %.2f)\n', ...
  Pa.E/1e9, Pa.E_lo/1e9, Pa.E_hi/1e9, E_true/1e9);
assert(abs(log(Pa.E/E_true)) < log(1.25), ...
  'patch invented structure on constant truth (%.2f vs %.2f GPa)', ...
  Pa.E/1e9, E_true/1e9);

% (b) patch on the flat shelf: moment and curvature both ~0, so the data
% cannot know the stiffness there. Unconstrained is the correct answer.
popts.E_patch = struct('lo', 11000, 'hi', 13000, 'E_ref', Rj.E);
Pb = vdef.invertElasticModulus(xo, wo, popts);
span_b = log10(Pb.E_hi/Pb.E_lo);
fprintf('12b. flat-shelf patch: [%.2f, %.2f] GPa, span %.2f decades\n', ...
  Pb.E_lo/1e9, Pb.E_hi/1e9, span_b);
assert(~isfinite(span_b) || span_b > 0.8, ...
  ['a patch on the flat shelf reported a %.2f-decade interval; nothing ' ...
   'constrains stiffness there and the fit should say so'], span_b);

% (c) a genuinely stiffened patch is recovered. Truth built through the
% same D-equivalence the fit uses, on the independent fine grid.
STIFF = 2.2;
hv = repmat(h_true, size(xt));
hv(xt >= 3000 & xt <= 5000) = h_true * STIFF^(1/3);
Ftp = vdef.beamFlexure(xt, hv, E_true, struct('A0', 1.0));
Wc = zeros(size(xo)); inc = xo >= x0_true;
Wc(inc) = interp1(Ftp.x, Ftp.w, xo(inc), 'linear');
wc = Wc + sig*randn(size(Wc));
dc = zeros(size(xstr));
for q = 1:numel(xstr)
  xsmp = xstr(q) + BW*uu;
  ws = zeros(size(xsmp)); inb = xsmp >= x0_true;
  ws(inb) = interp1(Ftp.x, Ftp.d2wdx2, xsmp(inb), 'linear');
  dc(q) = K2*mean(ws);
end
dc = dc + sig_s*randn(size(dc));
popts.strain.y = dc;                       % the stiffened truth's strain
popts.E_patch = struct('lo', 3000, 'hi', 5000, 'E_ref', E_true);
Pc = vdef.invertElasticModulus(xo, wc, popts);
fprintf('12c. stiffened patch: E* = %.2f GPa (truth %.2f = %.1fx background)\n', ...
  Pc.E/1e9, STIFF*E_true/1e9, STIFF);
assert(abs(log(Pc.E/(STIFF*E_true))) < log(1.3), ...
  'stiffened patch recovered as %.2f GPa against a truth of %.2f', ...
  Pc.E/1e9, STIFF*E_true/1e9);

%% Figure
fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

hf = figure('Visible','off','Position',[100 100 1000 380]);

subplot(1,3,1);
% The analytic curve goes down FIRST and thick, so the thin numeric one
% drawn over it leaves its dashes visible. Plotted the other way round the
% two agree to 2e-4 of A0 and the check reads as a single line.
plot(x/1e3, wan, 'k--', 'LineWidth', 3); hold on;
plot(x/1e3, F.w, '-', 'LineWidth', 1.2);
plot(x/1e3, Fv.w, '-', 'LineWidth', 1.5);
grid on; xlabel('Seaward distance (km)'); ylabel('Deflection (m)');
title('Forward model');
legend('analytic','uniform h','h 380-220 m','Location','SouthEast');

subplot(1,3,2);
plot(xo/1e3, wo, 'o', 'MarkerSize', 3); hold on;
plot(R.x_grid/1e3, R.w_grid, 'k-', 'LineWidth', 1.5);
grid on; xlabel('Seaward distance (km)'); ylabel('Deflection (m)');
title(sprintf('Recovered E* = %.1f GPa', R.E/1e9));
legend('data','fit','Location','SouthEast');

subplot(1,3,3);
% Misfit relative to its own minimum, in log10, so the trade-off valley is
% visible rather than one dark pixel at the optimum.
Jrel = log10(R.J / min(R.J(:)));
imagesc(log10(R.E_grid/1e9), R.x0_grid/1e3, Jrel);
set(gca,'YDir','normal'); hold on;
plot(log10(R.E/1e9), R.x0/1e3, 'w+', 'MarkerSize', 10, 'LineWidth', 2);
caxis([0 1]); colorbar;
xlabel('log_{10} E* (GPa)'); ylabel('Clamp position (km)');
title('Misfit surface');

out_fn = fullfile(fig_dir,'beam_flexure_test.png');
print(hf, out_fn, '-dpng', '-r120');
fprintf('Wrote %s\n', out_fn);

fprintf('\ntest_beam_flexure: all checks passed\n');
