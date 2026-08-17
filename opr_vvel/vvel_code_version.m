function v = vvel_code_version()
%VVEL_CODE_VERSION The deployed code version, for provenance stamping.
%   Returns the contents of code_version.txt at the project root - written
%   by opr_vvel/server/deploy.sh at rsync time, since the deployed tree is
%   not a git checkout and cannot answer 'git describe' itself. In a local
%   git checkout (development, CI) the file is absent and git is asked
%   directly. When neither works the stamp says so explicitly rather than
%   guessing: a product marked 'unknown' is a product built outside the
%   deploy path, which is itself worth knowing.

proj = fileparts(fileparts(mfilename('fullpath')));
fn = fullfile(proj, 'code_version.txt');
if exist(fn, 'file')
  fid = fopen(fn, 'r');
  v = strtrim(fgetl(fid));
  fclose(fid);
  return;
end
[st, out] = system(sprintf( ...
  'cd ''%s'' && git describe --always --dirty --tags 2>/dev/null', proj));
if st == 0 && ~isempty(strtrim(out))
  v = ['git-' strtrim(out)];
else
  v = 'unknown (built outside deploy.sh and not a git checkout)';
end
end
