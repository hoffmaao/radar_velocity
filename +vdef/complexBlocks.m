function [I, W] = complexBlocks(igram, coh, ref_bin, starts, block, coh_thr, rows)
%COMPLEXBLOCKS Surface-referenced, coherence-weighted COMPLEX block means.
%   [I, W] = COMPLEXBLOCKS(igram, coh, ref_bin, starts, block, coh_thr, rows)
%   references every column of a multilooked interferogram to its own
%   surface bin and averages the result along track into blocks, keeping
%   the field COMPLEX. Nothing is unwrapped: this is the input the
%   coherent tidal stack (vdef.tidalStack) works on, and it is what lets
%   that estimator stand independent of the unwrap-then-difference chain.
%
%   Referencing is by phase only - each column is multiplied by the unit
%   conjugate of its reference-bin value - so the surface phase is removed
%   exactly and the amplitude structure of the column is untouched. A
%   column with no usable reference contributes nothing. Samples below
%   coh_thr carry zero weight, as in vdef.blockAverage.
%
%   igram, coh   Nt x Nx multilooked interferogram and coherence
%   ref_bin      1 x Nx surface reference bin per column (NaN = none)
%   starts       first column of each block; block = columns per block
%   coh_thr      coherence below which a sample is dropped
%   rows         fast-time rows to keep (the depth selection)
%
%   Returns I (nz x nb) the weighted complex mean and W (nz x nb) the
%   weight sum, so a caller can weight blocks by how much data made them.
%
%   See also vdef.tidalStack, vdef.blockAverage, vdef.multilook.

[Nt, Nx] = size(igram);
rows = rows(:); nz = numel(rows); nb = numel(starts);
w = coh; w(coh < coh_thr | ~isfinite(coh)) = 0;
ref_bin = ref_bin(:).';
have = isfinite(ref_bin);
ref_bin(have) = min(max(round(ref_bin(have)), 1), Nt);
igr = zeros(nz, Nx, 'like', igram);
for c = 1:Nx
  if ~have(c), w(:,c) = 0; continue; end
  r = igram(ref_bin(c), c);
  if ~isfinite(r) || abs(r) < eps, w(:,c) = 0; continue; end
  igr(:,c) = igram(rows, c) .* (conj(r)/abs(r));
end
igr(~isfinite(igr)) = 0;
wr = w(rows, :);
I = zeros(nz, nb); W = zeros(nz, nb);
for b = 1:nb
  cols = starts(b):min(starts(b)+block-1, Nx);
  W(:,b) = sum(wr(:,cols), 2);
  I(:,b) = sum(igr(:,cols) .* wr(:,cols), 2) ./ max(W(:,b), eps);
end
end
