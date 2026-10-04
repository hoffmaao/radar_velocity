function tf = pairAligned(o)
%PAIRALIGNED Whether a vvel pair output is fast-time aligned, hence usable.
%
%   tf = vdef.pairAligned(o) for a loaded vvel_task output o.
%
%   Outputs carry `aligned`: true when the product was built WITHOUT
%   z-motion compensation (the passes are aligned at the product level by
%   multipass coregistration, and coalignment is not applied - see
%   vdef.zmotionApplied), or when coalignment removed the compensation's
%   misalignment on a product built with it. Older outputs carry only
%   coalign_applied, which meant exactly that for the compensated builds
%   they come from. An output with neither predates coalignment and counts
%   as usable, as it always did.
%
%   coalign_applied keeps its literal meaning, "a coalignment shift was
%   applied", so it is false on every surface-coupled pair. Loaders must
%   gate on this function, not on that field, or they discard every pair of
%   a surface-coupled build.

if isfield(o, 'aligned') && ~isempty(o.aligned)
  tf = logical(o.aligned);
elseif isfield(o, 'coalign_applied') && ~isempty(o.coalign_applied)
  tf = logical(o.coalign_applied);
else
  tf = true;
end
end
