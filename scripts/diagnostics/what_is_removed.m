%WHAT_IS_REMOVED Is per-column coalignment removing ice signal or geometry?
%
%   Written to answer a direct challenge (17 Aug 2026): the scalar-era
%   results showed a strong tidal response; is per-column coalignment
%   FILTERING OUT a real signal rather than removing an artefact?
%
%   Result on that date, GL3 and GL1 networks (67 + 70 pairs):
%     removed-vs-GPS-prediction correlation: median r = 0.99 both
%     slope: median 1.00 (GL3), 0.93 (GL1)
%     removed along-track variation: 84-91 mm rms of apparent range
%     NOT predicted by GPS: 10-12 mm rms - and the coalign estimator's own
%       noise (~0.25 ns) is ~21 mm of apparent range, so the unpredicted
%       residual is FULLY accounted for by measurement noise, leaving no
%       room for an ice contribution.
%   What per-column coalignment removes is deterministic platform
%   geometry at tide scale, 50-100x larger than any possible strain
%   signal, and predictable from GPS with no radar input at all.
%
% For every pair of the GL3 and GL1 networks, compare the ALONG-TRACK
% PROFILE of the shift coalignment removed (dtau_bulk_win, measured from
% the radar data) against the profile PREDICTED from GPS platform heights
% alone: -(ref_z_sec(x) - ref_z_ref(x))/(c/2), averaged over the same
% windows. The fit below leaves the slope FREE, so it absorbs whatever
% fraction of the erroneous compensation survives calibration; the
% measured slope is the quantity of record here. (artefact_vs_signal.m
% instead applies a fixed ALPHA = 1.03 surviving fraction by convention,
% which is why its predictor carries the factor and this one does not.)
% GPS knows nothing about ice strain. So:
%   - the GPS-predicted part of the removed field is platform geometry,
%     and removing it cannot cost any ice signal;
%   - only the residual NOT predicted by GPS could possibly contain ice.
% Report, per pair: correlation and slope of removed-vs-predicted, and
% the rms of the unpredicted residual in mm of apparent column change -
% the UPPER BOUND on signal that per-column coalignment could be eating.

mp_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass';
root   = '/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground';
c = 299792458; n_ice = 1.78;

for pnc = {'EAGER_2022_GL3','EAGER_2022_GL1'}
  pn = pnc{1};
  vdir = fullfile(root,'CSARP_vvel_net');
  L = load(fullfile(mp_dir, sprintf('%s_multipass03.mat', pn)), 'pass');
  f = dir(fullfile(vdir, [pn '_vvel_*.mat']));
  slopes = []; rs = []; resid_mm = []; removed_mm = [];
  for q = 1:numel(f)
    tok = regexp(f(q).name, ['^' regexptranslate('escape',pn) '_vvel_(\d+)_(\d+)\.mat$'], ...
      'tokens','once');
    if isempty(tok), continue; end
    o = load(fullfile(vdir, f(q).name), 'dtau_bulk_win','coalign_x_win', ...
      'coalign_applied','pass_idx_ref','pass_idx_sec');
    if ~o.coalign_applied || numel(o.dtau_bulk_win) < 6, continue; end
    i1 = o.pass_idx_ref; i2 = o.pass_idx_sec;
    zr = L.pass(i1).ref_z(:); zs = L.pass(i2).ref_z(:);
    pred = -(zs - zr)/(c/2);                    % [s] per column
    % average the prediction into the coalign windows
    xw = o.coalign_x_win; mw = o.dtau_bulk_win;
    pw = nan(size(xw));
    half = round(mean(diff(xw(~isnan(xw))))/2);
    if ~isfinite(half) || half < 1, half = 150; end
    for k = 1:numel(xw)
      if ~isfinite(xw(k)), continue; end
      lo = max(1, round(xw(k))-half); hi = min(numel(pred), round(xw(k))+half);
      pw(k) = mean(pred(lo:hi), 'omitnan');
    end
    ok = isfinite(mw) & isfinite(pw);
    if nnz(ok) < 6, continue; end
    % centre both: the along-track VARIATION is what per-column adds over
    % scalar, so the mean is removed from each before comparing
    a = mw(ok) - mean(mw(ok)); b = pw(ok) - mean(pw(ok));
    if std(b) == 0, continue; end
    p1 = polyfit(b, a, 1);
    m = corrcoef(a, b);
    slopes(end+1) = p1(1); rs(end+1) = m(1,2); %#ok<AGROW>
    res = a - p1(1)*b;
    % convert residual time to mm of apparent one-way range in ice
    resid_mm(end+1) = 1e3 * std(res) * c/(2*n_ice); %#ok<AGROW>
    removed_mm(end+1) = 1e3 * std(a) * c/(2*n_ice); %#ok<AGROW>
  end
  fprintf('\n===== %s: removed along-track shift vs GPS prediction (%d pairs) =====\n', ...
    pn, numel(rs));
  fprintf('correlation r:    median %.2f  (range %.2f .. %.2f)\n', ...
    median(rs), min(rs), max(rs));
  fprintf('slope:            median %.2f  (range %.2f .. %.2f)\n', ...
    median(slopes), min(slopes), max(slopes));
  fprintf('removed variation: median %.1f mm rms of apparent range\n', median(removed_mm));
  fprintf('NOT predicted by GPS: median %.1f mm rms  <- ceiling on any ice signal removed\n', ...
    median(resid_mm));
  clear L;
end
fprintf(['\nReading: r near 1 with slope near 1 means the removed field is\n' ...
  'deterministic platform geometry. The unpredicted residual is the most\n' ...
  'ice signal that could be hiding in it - compare with the ~0.4 mm/m x\n' ...
  'tide-range = ~0.4 mm scale of real strain over the coalign window, and\n' ...
  'with the 1.5 mm per-block noise.\n']);
