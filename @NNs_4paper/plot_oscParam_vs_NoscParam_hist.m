function plot_oscParam_vs_NoscParam_hist(obj,paramName,norm_flag)
% plot the parameter distribution if the oscillatory CPGs and the
% non-oscillatory CPGs.

% *) 'norm_flag' - '1' to norm parameter based on min/Max values

%number of bins to show:
binsNum = 20;

seq_osc = (vertcat(obj.results(obj.osc_ids).seq))';
% seq_osc = seq_osc(1:18,:);

seq_n_osc = (vertcat(obj.results(~obj.osc_ids).seq))';
% seq_n_osc = seq_n_osc(1:18,:);

switch paramName
    case {'all'}
        for i = 1:length(obj.seqOrder)
            p_name = obj.seqOrder{1,i};

            p_vec = seq_n_osc(strcmp(p_name,obj.seqOrder),:);
            p_osc_vec = seq_osc(strcmp(p_name,obj.seqOrder),:);
            
            if norm_flag
                p_vec = obj.norm_min_max(p_vec,p_name);
                p_osc_vec = obj.norm_min_max(p_osc_vec,p_name);
            end
            
            hist_compare(p_vec,p_osc_vec,p_name,binsNum,...
                {'n-osc CPGs','osc CPGs'},'plot');
        end
    otherwise
        p_vec = seq_n_osc(strcmp(paramName,obj.seqOrder),:);
        p_osc_vec = seq_osc(strcmp(paramName,obj.seqOrder),:);
        
        if norm_flag
            p_vec = obj.norm_min_max(p_vec,paramName);
            p_osc_vec = obj.norm_min_max(p_osc_vec,paramName);
            paramName = [paramName,' norm'];
        end
        
        hist_compare(p_vec,p_osc_vec,paramName,binsNum,...
                {'n-osc CPGs','osc CPGs'},'plot');
end

end

