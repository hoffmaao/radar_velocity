function d = figureDir(name)
%FIGUREDIR The one directory every figure and fitted product is written to.
%   d = FIGUREDIR() returns <repo root>/figs, creating it if it does not
%   exist. d = FIGUREDIR(name) returns the full path of a file in it.
%
%   WHY THIS EXISTS. Every figure script used to carry its own absolute
%   path to a scratch directory on the processing server, and there were
%   two of them - vvel/figures for the survey figures and
%   vvel/figures_flexure for the flexure ones - so a figure's home
%   depended on which script drew it, the outputs were split across two
%   directories that no checkout contained, and moving the analysis to
%   another machine meant editing a dozen scripts. The location is now
%   derived from this file's own position, so the figures land beside the
%   code that made them wherever that code sits: in the working copy when
%   run locally, in the deployed copy when run on the server. figs/ is
%   gitignored, so nothing here is ever committed.
%
%   The fitted products (flexure_fit_cats.mat, tidal_stack.mat,
%   tidal_depth_profile.mat) live here too rather than in a directory of
%   their own. They are the exact inputs the figures were drawn from, and
%   keeping them together is what lets a figure script redraw without
%   refitting: strain_flexure and tidal_stack both look for the flexure
%   fit in this directory by default.
%
%   Set the environment variable RADAR_VELOCITY_FIGS to send everything
%   somewhere else without touching a script - useful for a one-off run
%   whose output should not overwrite the current figure set.
%
%   See also scripts/figures/grl_figure.m.

d = getenv('RADAR_VELOCITY_FIGS');
if isempty(d)
  d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figs');
end
if ~exist(d, 'dir')
  [ok, msg] = mkdir(d);
  assert(ok, 'could not create the figure directory %s: %s', d, msg);
end
if nargin >= 1 && ~isempty(name)
  d = fullfile(d, name);
end
end
