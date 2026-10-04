function h_fig = get_figures(num_figs, visible_flag, user_data)
% Stub of OPR get_figures for standalone testing (always invisible)
h_fig = [];
for k = 1:num_figs
  hf = figure('Visible','off');
  if exist('theme', 'file'), theme(hf, 'light'); end   % not the OS dark mode (see grl_figure)
  h_fig(end+1) = hf; %#ok<AGROW>
end
end
