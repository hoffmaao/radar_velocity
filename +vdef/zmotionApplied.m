function tf = zmotionApplied(param_multipass)
%ZMOTIONAPPLIED Whether multipass motion-compensated each pass by ref_z.
%
%   tf = vdef.zmotionApplied(param_multipass) reads the setting a multipass
%   product was built with. multipass.m shifts every pass by ref_z/(c/2)
%   unless param.multipass.zmotion_comp_en is false, and products built
%   before that switch existed always applied it, so an absent field (or an
%   absent param_multipass) means true.
%
%   It matters downstream because the coalignment PREDICTION is built from
%   the compensation: with it applied, a secondary pass trails the
%   reference by -(ref_z_sec - ref_z_ref)/(c/2), which on a floating shelf
%   is the tidal heave mis-applied as range. Built without it, the antenna
%   rides the surface, the predicted misalignment is zero, and a measured
%   coalignment shift near zero is the expected outcome rather than a
%   failure.

tf = true;
if isempty(param_multipass) || ~isstruct(param_multipass), return; end
if isfield(param_multipass, 'multipass'), param_multipass = param_multipass.multipass; end
if isfield(param_multipass, 'zmotion_comp_en') && ~isempty(param_multipass.zmotion_comp_en)
  tf = logical(param_multipass.zmotion_comp_en);
end
end
