function [tide, pass_ok, info] = pass_tide(pn, mp_dir, opts)
%PASS_TIDE The tide regressor and pass gate for one multipass product.
%   [tide, pass_ok, info] = PASS_TIDE(pn, mp_dir, opts) returns, for every
%   pass of product pn, the CATS2008 tide at the pass mid-time (NaN for a
%   pass the gate rejected or that has no heights), and the gate itself.
%   It is the one place the tidal chain gets its regressor from, so the
%   surface admittance, the englacial strain admittance, and every figure
%   that regresses on the tide use the SAME tide and drop the SAME passes.
%
%   WHY NOT THE LINE-MEAN HEIGHT. The chain first used each pass's mean
%   GPS height as its tide. Against CATS2008 those means carry a 13-15 cm
%   rms non-tidal residual, half of it uniform along the line, and a
%   regressor carrying the same error as the data pulls every slope
%   toward one: the surface admittance profile was compressed, the strain
%   admittance attenuated by 15-25%, and a single pass with a 1.1 m height
%   error (GL2, 20221207_03) sat in every fit undetected. The mechanism
%   and its size are in vdef.surfaceAdmittance; the gate lives there too,
%   and this function only runs it on the product's ref_z.
%
%   opts (all optional): .block (200), .min_pass (5), .max_pass_sigma (4)
%                        - passed to vdef.surfaceAdmittance; .quad_period_h
%                        (24.84) for info.tide_q
%
%   info: .tmid   pass mid-times [UTC s since 1970]
%         .resid  line-mean residual about its fit to the tide [m]
%         .abar   slope of the line mean on the tide
%         .segs   day_seg of each pass, where the product records it
%         .main_idx
%         .tide_q  the QUADRATURE tide [m]: T/(2*pi) times the CATS2008 tide
%                  rate, T = opts.quad_period_h (24.84 h, the McMurdo diurnal
%                  band), NaN where tide is. For a diurnal tide it is the
%                  tide a quarter period on; vdef.tidalStack (opts.dq) uses
%                  it to test whether the response lags the tide
%
%   See also cats2008_tide, vdef.surfaceAdmittance,
%   scripts/diagnostics/elastic_modulus.m.

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'block')          || isempty(opts.block),          opts.block = 200;        end
if ~isfield(opts,'min_pass')       || isempty(opts.min_pass),       opts.min_pass = 5;       end
if ~isfield(opts,'max_pass_sigma') || isempty(opts.max_pass_sigma), opts.max_pass_sigma = 4; end
if ~isfield(opts,'quad_period_h')  || isempty(opts.quad_period_h),  opts.quad_period_h = 24.84; end
addpath(fileparts(mfilename('fullpath')));            % cats2008_tide

D = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass', 'param_multipass');
main_idx = D.param_multipass.multipass.baseline_master_idx;
Np = numel(D.pass); Nx = numel(D.pass(main_idx).ref_z);
Z = nan(Nx, Np); tmid = nan(1, Np); segs = repmat({''}, 1, Np);
for k = 1:Np
  z = D.pass(k).ref_z(:);
  if numel(z) ~= Nx, continue; end
  Z(:,k)  = z;
  tmid(k) = mean(D.pass(k).gps_time, 'omitnan');
  if isfield(D.pass(k),'param_pass') && isfield(D.pass(k).param_pass,'day_seg')
    segs{k} = char(D.pass(k).param_pass.day_seg);
  end
end
clear D;

cats = cats2008_tide(tmid);
S = vdef.surfaceAdmittance(Z, cats, struct('block', opts.block, ...
      'min_pass', opts.min_pass, 'max_pass_sigma', opts.max_pass_sigma));
pass_ok = S.pass_ok;
tide = cats; tide(~pass_ok) = NaN;
hs = 1800;                                            % s, central difference on the 5-min prediction
rate = (cats2008_tide(tmid + hs) - cats2008_tide(tmid - hs)) / (2*hs);
tide_q = opts.quad_period_h*3600/(2*pi) * rate; tide_q(~pass_ok) = NaN;
info = struct('tmid', tmid, 'resid', S.pass_resid, 'abar', S.abar, ...
              'segs', {segs}, 'main_idx', main_idx, 'tide_q', tide_q);
end
