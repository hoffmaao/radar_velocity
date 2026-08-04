function fn = opr_filename_out(param, out_path, tmp_dir, generic_data_flag)
% Stub of OPR opr_filename_out for standalone testing: paths under
% param.stub_out_root instead of gRadar.out_path.
%
% generic_data_flag excludes the day_seg directory, which is the mode both
% multipass and vvel use - a repeat-pass product spans segments, so it is
% filed by pass_name at the season level rather than under any one of them.
if nargin < 4 || isempty(generic_data_flag)
  generic_data_flag = 0;
end
fn = fullfile(param.stub_out_root, ['CSARP_' out_path]);
if ~generic_data_flag
  fn = fullfile(fn, param.day_seg);
end
end
