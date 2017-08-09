function [ output_args ] = plot_oscParam_vs_NoscParam( input_args )
% plot the parameter distribution if the oscillatory CPGs and the
% non-oscillatory CPGs. also, calc the "Kullback-Leibler" divergence

seq_osc = (vertcat(obj.results(obj.ids).seq))';
seq_osc = seq_osc(1:18,:);

seq_n_osc = (vertcat(results(:).seq))';
seq_n_osc = seq_n_osc(1:18,:);

for i = 1:length(obj.seqOrder)
    p_name = obj.seqOrder{1,i};
    
    p_vec = seq_n_osc(strcmp(p_name,obj.seqOrder),:);
    p_osc_vec = seq_osc(strcmp(p_name,obj.seqOrder),:);
    
    hist_compare(p_vec,p_osc_vec,p_name,20,{'n-osc CPGs','osc CPGs'},'plot');

end

clear p_name p_vec p_osc_vec pvalue rejection

% % % % % "Kullback-Leibler" divergence:
clc

tau_n_osc = seq_n_osc(strcmp('tau',obj.seqOrder),:);
tau_osc = seq_osc(strcmp('tau',obj.seqOrder),:);
dist_tau = KL_div_4paper(tau_n_osc,tau_osc)
clear tau_n_osc tau_osc

b_n_osc = seq_n_osc(strcmp('b',obj.seqOrder),:);
b_osc = seq_osc(strcmp('b',obj.seqOrder),:);
dist_b = KL_div_4paper(b_n_osc,b_osc)
clear b_n_osc b_osc


end

