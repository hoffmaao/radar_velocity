function AP = read_apres_xy(fn)
%READ_APRES_XY Site name and EPSG:3031 metres from the whitespace table
%   `x  y  name`. Returns an empty struct (with a note) when the file is
%   absent, so the figure still builds without the site coordinates.
AP = struct('name',{},'x',{},'y',{},'along',{},'off',{});
if ~exist(fn,'file')
  fprintf('ApRES site coordinates not found (%s) - sites not drawn.\n', fn);
  return;
end
try
  fid = fopen(fn,'r');
  C = textscan(fid, '%f %f %s');
  fclose(fid);
  for k = 1:numel(C{3})
    AP(end+1) = struct('name', C{3}{k}, 'x', C{1}(k), 'y', C{2}(k), ...
      'along', NaN, 'off', NaN); %#ok<AGROW>
  end
  fprintf('ApRES sites read: %s\n', strjoin({AP.name}, ', '));
catch ME
  fprintf('Could not read %s (%s) - sites not drawn.\n', fn, ME.message);
end
end
