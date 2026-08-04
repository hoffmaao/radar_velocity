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
%   Returns:
%     blk.dtau      Nt x Nblk coherence-weighted mean traveltime difference
%     blk.dtau_std  Nt x Nblk weighted standard deviation within the block
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

  sw = sum(wb, 2);
  cover = sum(wb > 0, 2) / n;
  blk.coverage(:,b) = cover;

  keep = sw > 0 & cover >= opts.min_coverage;

  mu = sum(wb .* db, 2) ./ max(sw, eps);
  blk.dtau(keep, b) = mu(keep);

  % Weighted variance about the weighted mean
  var_num = sum(wb .* bsxfun(@minus, db, mu).^2, 2);
  sd = sqrt(var_num ./ max(sw, eps));
  blk.dtau_std(keep, b) = sd(keep);

  cb = coh(:,cols);
  cb(wb == 0) = NaN;
  mc = mean(cb, 2, 'omitnan');
  blk.coh(keep, b) = mc(keep);

  blk.Surface(b) = mean(map.Surface(cols), 'omitnan');
end

end
