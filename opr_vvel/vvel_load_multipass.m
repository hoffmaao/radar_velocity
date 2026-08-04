function mp = vvel_load_multipass(param)
% mp = vvel_load_multipass(param)
%
% Loads a multipass comp_mode 3 product and resolves the two mappings the
% rest of the module needs but the file does not state directly:
%
%   mp.baseline_main_idx  the pass every image was coregistered onto, so
%                           its time / along_track / surface vectors are
%                           the common axes
%   mp.pass_en_idxs         data(:,:,k) holds pass(pass_en_idxs(k)) -
%                           multipass concatenates only the ENABLED passes
%                           into data while leaving pass complete, so the
%                           two index spaces differ whenever pass_en_mask
%                           has a false in it
%
% Returns [] (with a warning) when the file is missing or incomplete, so
% callers can skip rather than abort a batch.
%
% mp fields: data, pass, in_fn, baseline_main_idx, pass_en_idxs, and
% param_multipass / param_combine_passes when present.
%
% See also: vvel.m, vvel_task.m, multipass.m

mp = [];

%% Resolve the input file
% =========================================================================
if isfield(param.vvel,'in_fn') && ~isempty(param.vvel.in_fn)
  % Absolute path override: useful when running against a product copied
  % off the server without the surrounding CSARP_ directory structure
  in_fn = param.vvel.in_fn;
else
  in_dir = opr_filename_out(param,param.vvel.in_path,'',1);
  in_fn = fullfile(in_dir, sprintf('%s%s_multipass%02.0f.mat', ...
    param.vvel.pass_name, param.vvel.output_fn_midfix, param.vvel.comp_mode));
end

if ~exist(in_fn,'file')
  warning('The multipass file does not exist. Perhaps param.vvel.in_path or param.vvel.pass_name is wrong, or multipass.m has not been run with comp_mode 3. File does not exist:\n  %s.', in_fn);
  return;
end

%% Load
% =========================================================================
% Only data and pass are read. ref is a copy of the main pass taken
% before the coregistration loop and can still carry a full SLC image
% (GBs), and everything it holds that this module uses is also in
% pass(baseline_main_idx).
have = whos('-file', in_fn);
required = {'data','pass'};
missing = setdiff(required, {have.name});
if ~isempty(missing)
  warning('Required variable(s) %s missing from the multipass file. Perhaps it was written with a comp_mode other than 3. File:\n  %s.', ...
    strjoin(missing, ', '), in_fn);
  return;
end

want = {'data','pass','param_multipass','param_combine_passes'};
sel  = intersect(want, {have.name});
fprintf('Loading multipass product:\n  %s\n', in_fn);
mp = load(in_fn, sel{:});
mp.in_fn = in_fn;

if ndims(mp.data) < 3
  warning('The multipass product holds a single pass image; a repeat-pass pair needs at least two. File:\n  %s.', in_fn);
  mp = [];
  return;
end

%% Baseline main pass and the data-slice mapping
% =========================================================================
Npass = numel(mp.pass);
Nslice = size(mp.data,3);

baseline_main_idx = 1;   % multipass default
pass_en_mask = [];
if isfield(mp,'param_multipass') && isfield(mp.param_multipass,'multipass')
  pm = mp.param_multipass.multipass;
  if isfield(pm,'baseline_master_idx') && ~isempty(pm.baseline_master_idx)
    baseline_main_idx = pm.baseline_master_idx;
  end
  if isfield(pm,'pass_en_mask') && ~isempty(pm.pass_en_mask)
    pass_en_mask = logical(pm.pass_en_mask);
  end
end
if ~isempty(param.vvel.baseline_main_idx)
  baseline_main_idx = param.vvel.baseline_main_idx;
end

if isempty(pass_en_mask)
  pass_en_idxs = 1:Npass;
else
  % multipass treats any pass beyond the end of the mask as enabled
  pass_en_mask(end+1:Npass) = true;
  pass_en_idxs = find(pass_en_mask(1:Npass));
end

if numel(pass_en_idxs) ~= Nslice
  error('vvel:pass_en_mask', ...
    ['The enabled-pass mask names %d passes but data holds %d images, so slices cannot be attributed to passes. ' ...
     'Check param_multipass.multipass.pass_en_mask in:\n  %s.'], ...
    numel(pass_en_idxs), Nslice, in_fn);
end

if baseline_main_idx < 1 || baseline_main_idx > Npass
  error('vvel:baseline_main_idx','baseline_main_idx %d is outside the %d passes in:\n  %s.', ...
    baseline_main_idx, Npass, in_fn);
end

% The main pass supplies the common axes AND is one half of every pair the
% default 'main' pairing builds, so a disabled main is a dead batch. Say so
% once here rather than once per pair in vvel_task.
if ~ismember(baseline_main_idx, pass_en_idxs)
  error('vvel:baseline_main_idx', ...
    ['baseline_main_idx %d is disabled in pass_en_mask, so it has no image in data ' ...
     'and cannot be one half of a pair. Enable it, or point param.vvel.baseline_main_idx ' ...
     'at an enabled pass, in:\n  %s.'], baseline_main_idx, in_fn);
end

mp.baseline_main_idx = baseline_main_idx;
mp.pass_en_idxs = pass_en_idxs;

fprintf('  %d passes (%d enabled), baseline main %d, %d x %d x %d image stack\n', ...
  Npass, Nslice, baseline_main_idx, size(mp.data,1), size(mp.data,2), Nslice);

end
