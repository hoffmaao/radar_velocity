function [h, S] = grl_figure(width_mm, height_mm)
%GRL_FIGURE Create a figure sized and fonted for an AGU/GRL submission.
%
%   [h, S] = GRL_FIGURE(width_mm, height_mm) returns an invisible white
%   figure whose PHYSICAL size is width_mm x height_mm and whose default
%   fonts are S.font at S.base_pt points, plus a struct S of the style
%   constants for anything that has to be set explicitly.
%
%   WHY PHYSICAL SIZE. AGU asks for lettering close to 8 point "at final
%   print size" (6 point for sub/superscript), in Arial, Helvetica, Times
%   or Symbol, in a figure 50-170 mm wide - one column 50-85 mm, two
%   column 105-170 mm. "At final print size" is the whole difficulty: a
%   figure authored 1150 px wide and printed at -r140 is 209 mm, so
%   dropping it into a 170 mm column scales every font by 0.81, and a
%   figure authored at a different pixel size scales by something else.
%   Sizes set in the script then have no fixed relationship to the page
%   and no two figures in the set match. Authoring at print size removes
%   the scaling step: what the script sets is what the page gets.
%
%   Sets PaperPosition to match, so `print -dpng -r<N>` rasterises the
%   same physical figure at N dots per inch - resolution changes, size
%   does not.
%
%   Axes positions in these scripts are NORMALISED, so a figure resized
%   this way keeps its layout; only the physical scale and the fonts
%   change.
%
%   Usage:
%     [h, S] = grl_figure(S_TWO_COL, 105);
%     ax = axes('parent',h,'Position',[0.1 0.1 0.8 0.8]);
%     xlabel(ax,'Along track (km)');            % inherits 8 pt Helvetica
%     print(h, fn, '-dpng', sprintf('-r%d', S.dpi));

S = grl_constants();
if nargin < 1 || isempty(width_mm),  width_mm  = S.two_col_mm; end
if nargin < 2 || isempty(height_mm), height_mm = width_mm*0.62; end

if width_mm < S.min_mm || width_mm > S.max_mm
  error('grl_figure:width', ...
    'width %.1f mm is outside AGU''s %.0f-%.0f mm range', ...
    width_mm, S.min_mm, S.max_mm);
end

wcm = width_mm/10; hcm = height_mm/10;
h = figure('Visible','off', 'Color','w', ...
  'Units','centimeters', 'Position',[2 2 wcm hcm], ...
  'PaperUnits','centimeters', 'PaperPosition',[0 0 wcm hcm], ...
  'PaperSize',[wcm hcm]);
% Pin the light theme. From R2025a a figure follows the OPERATING SYSTEM's
% appearance, so on a Mac in dark mode (automatic at night) every axes
% exports with a black background and grey text, while explicitly coloured
% elements stay as set - a half-dark figure. 'Color','w' above does not
% prevent it. R2024b (the server) has no themes and needs nothing.
if exist('theme', 'file'), theme(h, 'light'); end

% Defaults must be set on the FIGURE before any axes are created; setting
% them afterwards leaves the existing axes on MATLAB's own 10 pt default
% and produces a figure that is half compliant.
set(h, ...
  'DefaultAxesFontName',      S.font, ...
  'DefaultAxesFontSize',      S.base_pt, ...
  'DefaultTextFontName',      S.font, ...
  'DefaultTextFontSize',      S.base_pt, ...
  'DefaultLegendFontName',    S.font, ...
  'DefaultLegendFontSize',    S.base_pt, ...
  'DefaultColorbarFontName',  S.font, ...
  'DefaultColorbarFontSize',  S.base_pt, ...
  'DefaultAxesLineWidth',     S.axis_lw, ...
  'DefaultLineLineWidth',     S.line_lw);
end

%% ========================================================================
function S = grl_constants()
%GRL_CONSTANTS The numbers, in one place so the .m and .py agree.
S.one_col_mm = 85;      % AGU one-column figure: 50-85 mm
S.two_col_mm = 170;     % AGU two-column figure: 105-170 mm
S.min_mm     = 50;
S.max_mm     = 170;
S.base_pt    = 8;       % lettering close to 8 pt at final print size
S.small_pt   = 7;       % dense in-plot annotation only, never axis text
S.sub_pt     = 6;       % AGU's floor, for sub/superscript
% MATLAB maps 'Helvetica' to a Helvetica-metric sans on every platform we
% build on; it is on AGU's accepted list, as is Arial.
S.font       = 'Helvetica';
S.axis_lw    = 0.6;
S.line_lw    = 1.2;
S.dpi        = 400;     % raster output; type is rasterised, so go high
end
