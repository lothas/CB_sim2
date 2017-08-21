function plot_2D_hist(obj,paramNames,norm_flag)
% plot the parameter distribution if the oscillatory CPGs and the
% non-oscillatory CPGs.

% *) 'paramNames' - 1X2 cell array with the names of the params
%number of bins to show:
binsNum = 20;

% oscillating CPGs:
seq_osc = (vertcat(obj.results(obj.osc_ids).seq))';
arr = seq_osc(1:18,:);

% % oscillating CPGs (in range):
% seq_osc_inRange = (vertcat(obj.results(obj.osc_inRange_ids).seq))';
% arr = seq_osc_inRange(1:18,:);

% %  non oscillating CPGs:
% seq_n_osc = (vertcat(obj.results(~obj.osc_ids).seq))';
% arr = seq_n_osc(1:18,:);

vec1 = arr(strcmp(paramNames{1,1},obj.seqOrder),:);
vec2 = arr(strcmp(paramNames{1,2},obj.seqOrder),:);

if norm_flag
    vec1 = obj.norm_min_max(vec1,paramNames{1,1});
    vec2 = obj.norm_min_max(vec2,paramNames{1,2});
    paramNames{1,1} = [paramNames{1,1},' norm'];
    paramNames{1,2} = [paramNames{1,2},' norm'];
end

histogram2(vec1,vec2,binsNum,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of \tau and b');
xlabel(paramNames{1,1});
ylabel(paramNames{1,2});

end