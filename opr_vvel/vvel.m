function [success] = vvel(param,param_override)
% [success] = vvel(param,param_override)
%
% Infers the vertical strain rate eps_zz versus depth and along-track
% position from a repeat-pass CSARP_multipass comp_mode 3 product
% (multipass.m output). For each requested pass PAIR the interferogram
% phase is converted into a differential two-way traveltime referenced to
% the surface return, averaged along track, converted into the vertical
% displacement of each reflector relative to the surface, and inverted for
% a Legendre strain-rate profile (see the +vdef package, which must be on
% the MATLAB path).
%
% Output .mat files (CSARP_vvel) and image files are written, one set per
% pair.
%
% WHY THIS IS NOT A CLUSTER JOB. Unlike the per-frame OPR products, a
% multipass product is a single file holding every pass in the experiment,
% and loading it is the dominant cost of the run. One cluster task per pair
% would reload the whole stack for every pair. So vvel.m loads it once and
% runs the pairs in-process, which also matches how multipass.m itself is
% invoked (directly, no cluster). vvel_task.m still takes a single pair and
% still returns a success flag, so it can be wrapped in a cluster task
% later if a season ever gets big enough to need one.
%
% See also: run_vvel.m, vvel_task.m, vvel_load_multipass.m, multipass.m, +vdef

%% General Setup
% =====================================================================
if exist('param_override','var')
  param = merge_structs(param, param_override);
end
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.vvel.pass_name, datestr(now));
fprintf('=====================================================================\n');

%% Input arguments check and setup
% =========================================================================

% er_ice, c = speed of light
physical_constants;
[output_dir,radar_type,radar_name] = opr_output_dir(param.radar_name);

param = vvel_defaults(param);

%% Load the product once, resolve the pairs
% =====================================================================
mp = vvel_load_multipass(param);
if isempty(mp)
  success = false;
  return;
end

pairs = vvel_resolve_pairs(param, mp);
if isempty(pairs)
  warning('No pass pairs to process for %s. Check param.vvel.pairs and the enabled-pass mask.', param.vvel.pass_name);
  success = false;
  return;
end
fprintf('%d pair(s) to process\n', size(pairs,1));

%% Process each pair
% =====================================================================
out_dir = opr_filename_out(param,param.vvel.out_path,'',1);
success = true;
for pair_idx = 1:size(pairs,1)
  ref_idx = pairs(pair_idx,1);
  sec_idx = pairs(pair_idx,2);

  fn_name = sprintf('%s%s_vvel_%02d_%02d', param.vvel.pass_name, ...
    param.vvel.output_fn_midfix, ref_idx, sec_idx);
  out_fn = fullfile(out_dir,[fn_name '.mat']);

  if param.vvel.rerun_only && exist(out_fn,'file')
    fprintf('  Already exists [rerun_only skipping]: %s\n', out_fn);
    continue;
  end

  fprintf('---------------------------------------------------------------------\n');
  fprintf('Pair %d of %d: pass %d -> %d (%s)\n', ...
    pair_idx, size(pairs,1), ref_idx, sec_idx, datestr(now));

  param_pair = param;
  param_pair.load.pair = [ref_idx sec_idx];

  if param.vvel.continue_on_error
    % A single bad pair (an incoherent baseline, a pass with no gps_time)
    % should not cost the rest of the batch
    try
      pair_success = vvel_task(param_pair, mp);
    catch ME
      warning('Pair %d -> %d failed: %s', ref_idx, sec_idx, ME.message);
      pair_success = false;
    end
  else
    pair_success = vvel_task(param_pair, mp);
  end
  success = success && pair_success;
end

fprintf('Done %s\n', datestr(now));

end

%% ========================================================================
function pairs = vvel_resolve_pairs(param, mp)
% Turn param.vvel.pairs into an N x 2 list of [ref sec] pass indices.
% Only enabled passes can be paired - a disabled pass has no image.

en = mp.pass_en_idxs(:).';
main_pass = mp.baseline_main_idx;

if isnumeric(param.vvel.pairs)
  pairs = param.vvel.pairs;
  if size(pairs,2) ~= 2
    error('vvel:pairs','param.vvel.pairs must be an N x 2 matrix of [ref sec] pass indices.');
  end
  bad = ~ismember(pairs(:,1), en) | ~ismember(pairs(:,2), en);
  if any(bad)
    warning('Dropping %d requested pair(s) that name a pass with no image in the product (disabled in pass_en_mask).', nnz(bad));
    pairs = pairs(~bad,:);
  end
else
  switch lower(param.vvel.pairs)
    case 'main'
      others = en(en ~= main_pass);
      pairs = [repmat(main_pass, numel(others), 1) others(:)];
    case 'sequential'
      pairs = [en(1:end-1).' en(2:end).'];
    case 'all'
      pairs = zeros(0,2);
      for a = 1:numel(en)
        for b = a+1:numel(en)
          pairs(end+1,:) = [en(a) en(b)]; %#ok<AGROW>
        end
      end
    otherwise
      error('vvel:pairs','Unrecognised param.vvel.pairs "%s"; use ''main'', ''sequential'', ''all'', or an N x 2 matrix.', param.vvel.pairs);
  end
end

self = pairs(:,1) == pairs(:,2);
if any(self)
  warning('Dropping %d self-pair(s).', nnz(self));
  pairs = pairs(~self,:);
end

end
