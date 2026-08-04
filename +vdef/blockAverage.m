function blk = blockAverage(dtau, map, info, opts)
%BLOCKAVERAGE Coherence-weighted along-track averaging of dtau.
%   blk = BLOCKAVERAGE(dtau, map, info, opts) averages the differential
%   traveltime map into along-track blocks of opts.block_size range lines,
%   weighting each sample by its coherence and rejecting bins whose
%   coherent coverage within the block falls below opts.min_coverage.
%
%   Averaging along track before inverting is what buys the precision: the
%   per-pixel phase noise at coherence gamma is roughly
%   sigma_phi ~ sqrt((1-gamma^2)/(2*gamma^2*L)) for L looks, so a block of
%   ~1000 range lines takes a 0.1 rad per-pixel scatter down to well under
%   a milliradian, i.e. sub-millimetre displacement at 750 MHz.
%
%   EFFECTIVE SAMPLE COUNT. The scatter of the samples inside a block is
%   not the uncertainty of the block mean - averaging n independent samples
%   divides it by sqrt(n) - and that is the whole point of averaging here.
%   The samples are not independent, though: vdef.multilook has already run
%   a boxcar of opts.mlook_window(2) range lines along track, so adjacent
%   columns share most of their input. The count is therefore deflated by
%   opts.cols_per_look (default: the along-track multilook length), the
%   along-track counterpart of the bins_per_look correction that
%   vdef.invertStrainRate applies in fast time. Unequal coherence weights
%   are handled with the Kish effective count, (sum w)^2 / sum(w^2).
%
%   Returns:
%     blk.dtau      Nt x Nblk coherence-weighted mean traveltime difference
%     blk.dtau_std  Nt x Nblk 1-sigma of that MEAN: the within-block
%                   weighted standard deviation divided by sqrt(n_eff)
%     blk.dtau_scatter  Nt x Nblk weighted standard deviation WITHIN the
%                   block; a relative noise measure, not an error bar
%     blk.n_eff     Nt x Nblk effective independent samples behind the mean
%     blk.coh       Nt x Nblk mean coherence of contributing samples
%     blk.coverage  Nt x Nblk fraction of the block that contributed
%     blk.cols      1 x Nblk cell array of column indices per block
%     blk.starts    1 x Nblk first column of each block
%     blk.Surface   1 x Nblk mean surface twtt per block
%
%   See also vdef.differentialRange, vdef.verticalDisplacement.

[Nt, Nx] = size(dtau);

block_size = max(1, round(opts.block_size));
starts = 1:block_size:Nx;
Nblk   = numel(starts);

% Along-track decorrelation length in columns
if isfield(opts,'cols_per_look') && ~isempty(opts.cols_per_look)
  cols_per_look = max(1, opts.cols_per_look);
elseif isfield(opts,'mlook_window') && numel(opts.mlook_window) >= 2
  cols_per_look = max(1, round(opts.mlook_window(2)));
else
  cols_per_look = 1;
end

coh = map.coherence;
if isempty(coh)
  coh = ones(Nt, Nx);
end
w = coh;
w(coh < opts.coherence_threshold) = 0;
w(~isfinite(dtau)) = 0;

d = dtau;
d(~isfinite(d)) = 0;

blk = [];
blk.dtau     = nan(Nt, Nblk);
blk.dtau_std = nan(Nt, Nblk);
blk.dtau_scatter = nan(Nt, Nblk);
blk.n_eff    = nan(Nt, Nblk);
blk.coh      = nan(Nt, Nblk);
blk.coverage = zeros(Nt, Nblk);
blk.cols     = cell(1, Nblk);
blk.starts   = starts;
blk.Surface  = nan(1, Nblk);

for b = 1:Nblk
  cols = starts(b):min(starts(b)+block_size-1, Nx);
  blk.cols{b} = cols;
  n = numel(cols);

  wb = w(:,cols);
  db = d(:,cols);

  sw  = sum(wb, 2);
  sw2 = sum(wb.^2, 2);
  cover = sum(wb > 0, 2) / n;
  blk.coverage(:,b) = cover;

  keep = sw > 0 & cover >= opts.min_coverage;

  mu = sum(wb .* db, 2) ./ max(sw, eps);
  blk.dtau(keep, b) = mu(keep);

  % Weighted variance about the weighted mean
  var_num = sum(wb .* bsxfun(@minus, db, mu).^2, 2);
  sd = sqrt(var_num ./ max(sw, eps));

  % Kish effective count of the contributing columns, deflated by the
  % along-track look length, then the standard error of the block mean
  n_eff = max((sw.^2 ./ max(sw2, eps)) / cols_per_look, 1);

  blk.dtau_scatter(keep, b) = sd(keep);
  blk.n_eff(keep, b)        = n_eff(keep);
  blk.dtau_std(keep, b)     = sd(keep) ./ sqrt(n_eff(keep));

  cb = coh(:,cols);
  cb(wb == 0) = NaN;
  mc = mean(cb, 2, 'omitnan');
  blk.coh(keep, b) = mc(keep);

  blk.Surface(b) = mean(map.Surface(cols), 'omitnan');
end

end
