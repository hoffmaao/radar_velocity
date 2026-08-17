%VERIFY_MULTIPASS_RERUN Prove a rerun multipass product matches the archive.
%
%   The reproducibility gate for the whole workflow: after any rerun of
%   opr_vvel/server/run_multipass_scratch.m, this compares the scratch
%   product against the archived one, variable by variable, and requires
%   the data cube to be BIT-IDENTICAL.
%
%   What "identical" has to tolerate, learned on the 4 Aug 2026 rerun of
%   GL1-GL4 (all four passed):
%     - pass gains fields the archive predates (toolbox 6fcec75f added
%       ref_roll/ref_pitch/ref_heading): NEW fields are reported, allowed
%     - proj_x/proj_y differ at ~1e-10 m (map reprojection rounding):
%       differences below GEO_TOL are allowed
%     - anything else differing is a FAILURE
%
%   Define PRODUCT before running; MASTER_OVERRIDE builds (suffixed _mNN)
%   have no archived counterpart and are refused rather than "verified"
%   against the wrong file.
%
%     matlab -batch "PRODUCT='EAGER_2022_GL3'; run('.../verify_multipass_rerun.m')"

if ~exist('PRODUCT','var')
  error('Define PRODUCT, e.g. PRODUCT=''EAGER_2022_GL3''');
end
PRODUCT = char(PRODUCT);
assert(isempty(regexp(PRODUCT,'_m\d\d$','once')), ...
  '%s is a master-override build; there is no archived counterpart to verify against', PRODUCT);
GEO_TOL = 1e-6;   % [m] projection rounding

arch = fullfile('/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_multipass', ...
  [PRODUCT '_multipass03.mat']);
scr  = fullfile('/kucresis/scratch/hoffmana_sta/vvel/2022_Antarctica_Ground/CSARP_multipass', ...
  [PRODUCT '_multipass03.mat']);
assert(exist(arch,'file')==2, 'no archive: %s', arch);
assert(exist(scr,'file')==2,  'no rerun: %s', scr);

A = load(arch); S = load(scr);
fail = false;

va = sort(fieldnames(A)); vs = sort(fieldnames(S));
if ~isequal(va, vs)
  fprintf('top-level variables differ:\n  archive-only: %s\n  rerun-only: %s\n', ...
    strjoin(setdiff(va,vs),','), strjoin(setdiff(vs,va),','));
  fail = true;
else
  fprintf('top-level vars: identical (%s)\n', strjoin(va.',', '));
end

% the data cube must be bit-identical
if ~isequaln(A.data, S.data)
  fprintf('FAIL: data cube differs\n'); fail = true;
else
  fprintf('data: BIT-IDENTICAL, size [%s]\n', num2str(size(S.data)));
end

% pass struct: same count required, new fields allowed, common fields must
% match (geo to tol). A rerun with a different pass count is the clearest
% possible non-reproduction (the EAGER_2022 10-vs-13 pathology), so it has
% to produce the clean verdict rather than an index error.
if numel(A.pass) ~= numel(S.pass)
  fprintf('FAIL: pass count differs: archive %d, rerun %d\n', ...
    numel(A.pass), numel(S.pass));
  fail = true;
else
  fa = fieldnames(A.pass); fs = fieldnames(S.pass);
  newf = setdiff(fs, fa); missf = setdiff(fa, fs);
  if ~isempty(newf),  fprintf('pass fields only in rerun (allowed): %s\n', strjoin(newf.',', ')); end
  if ~isempty(missf), fprintf('FAIL: pass fields LOST in rerun: %s\n', strjoin(missf.',', ')); fail = true; end
  bad = {};
  for c = intersect(fa, fs).'
    fn = c{1};
    for k = 1:numel(A.pass)
      x = A.pass(k).(fn); y = S.pass(k).(fn);
      if isequaln(x, y), continue; end
      if any(strcmp(fn, {'proj_x','proj_y'})) && isnumeric(x) && isequal(size(x),size(y)) ...
          && max(abs(x(:)-y(:))) < GEO_TOL
        continue;
      end
      bad{end+1} = sprintf('%s(pass %d)', fn, k); %#ok<AGROW>
      break;
    end
  end
  if isempty(bad)
    fprintf('pass: common fields match (geo fields to %.0e m)\n', GEO_TOL);
  else
    fprintf('FAIL: pass fields differ: %s\n', strjoin(bad,', ')); fail = true;
  end
end

if fail
  error('verify_multipass_rerun: %s DID NOT REPRODUCE', PRODUCT);
end
fprintf('\n%s: rerun REPRODUCES the archive.\n', PRODUCT);
