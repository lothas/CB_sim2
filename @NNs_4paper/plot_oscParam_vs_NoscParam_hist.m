function plot_oscParam_vs_NoscParam_hist(obj,paramName)
% plot the parameter distribution if the oscillatory CPGs and the
% non-oscillatory CPGs.

%number of bins to show:
binsNum = 20;

seq_osc = (vertcat(obj.results(obj.osc_ids).seq))';
seq_osc = seq_osc(1:18,:);

seq_n_osc = (vertcat(obj.results(~obj.osc_ids).seq))';
seq_n_osc = seq_n_osc(1:18,:);

switch paramName
    case {'all'}
        for i = 1:length(obj.seqOrder)
            p_name = obj.seqOrder{1,i};

            p_vec = seq_n_osc(strcmp(p_name,obj.seqOrder),:);
            p_osc_vec = seq_osc(strcmp(p_name,obj.seqOrder),:);

            hist_compare(p_vec,p_osc_vec,p_name,binsNum,...
                {'n-osc CPGs','osc CPGs'},'plot');
        end
    otherwise
        p_vec = seq_n_osc(strcmp(paramName,obj.seqOrder),:);
        p_osc_vec = seq_osc(strcmp(paramName,obj.seqOrder),:);
        hist_compare(p_vec,p_osc_vec,paramName,binsNum,...
                {'n-osc CPGs','osc CPGs'},'plot');
end

end

