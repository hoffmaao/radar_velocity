function [lines, dup] = surveyLines()
%SURVEYLINES The independent lines of the 2022 survey, and the duplicate.
%   lines = SURVEYLINES() returns the FOUR independent profiles.
%   [lines, dup] = SURVEYLINES() also returns the duplicate build, for the
%   few places that deliberately want it.
%
%   THERE ARE FOUR LINES, NOT FIVE. The multipass directory holds five
%   products, and EAGER_2022 was long carried alongside GL1 to GL4 as
%   though it were a fifth profile. It is not. It is a second build of
%   GL1, and the evidence is not subtle:
%     * the two centre lines are a median 0.5 m apart, against 140 m from
%       GL1 to its true neighbour GL2, 267 m to GL3 and 402 m to GL4;
%     * identical column count (2013), identical end points to five
%       decimal places, identical length (5.03 km);
%     * GL1's pass list is a strict superset - EAGER_2022's thirteen day
%       segments plus 20221209_01.
%   So the merge of the two IS GL1, which is why GL1 is the one to keep.
%   scripts/figures/leg1_merge_check.m reached this conclusion long ago
%   and said so in its header; the product lists were simply never
%   updated, so every table, map and pooled statistic counted leg 1 twice.
%
%   WHY THE DUPLICATE IS STILL WORTH HAVING. Two builds of the same ice
%   with almost the same passes bound the method's own systematic error in
%   a way no formal covariance can: their disagreement IS the error. That
%   is what scripts/figures/strain_rates.m and scripts/figures/tidal_evidence.m
%   use it for, what between_build.m measures, and what leg1_merge_check.m
%   tests the coalignment fix against. Those uses are legitimate and stay.
%   What is not legitimate is listing it beside GL1 to GL4 as though it
%   added an independent profile.
%
%   READ THE DISAGREEMENT BEFORE QUOTING EITHER. In the block-length sweep
%   the two builds of this one line do not merely differ in size: GL1 came
%   out consistent with zero at 100 m (-0.4 to +0.4 mm per m of tide,
%   sigma 0.6 to 1.2) while EAGER_2022 came out at -2.9 to -3.9 with the
%   same sign at every block length. Same ice, same track, differing by
%   one pass and by which pass is the reference. Whatever that number is
%   measuring, it is not only the ice.
%
%   See also vdef.surfaceAdmittance, scripts/figures/leg1_merge_check.m.

lines = {'EAGER_2022_GL1', 'EAGER_2022_GL2', 'EAGER_2022_GL3', 'EAGER_2022_GL4'};
dup   = struct('name', 'EAGER_2022', 'same_as', 'EAGER_2022_GL1', ...
               'note', 'second build of GL1, thirteen of its fourteen passes');
end
