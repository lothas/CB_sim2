function plot_osc_nosc_hist(subPlot_grid,seq,periods,seqOrder,osc_ids,n_osc_ids)
% PLOT_PARAM_HIST plot the histograms of all of the parameters
% 
% 'subPlot_grid' - 1x2 vector contain the subplot grid size =
% (N_param+period)
% 
N = length(seqOrder);

figure;
for i=1:(N+1)
    subplot(subPlot_grid(1,1),subPlot_grid(1,2),i);
    if (i > N)
        histogram(periods,100,'Normalization','pdf');
        title(['histogram of: ','periods']);
    else
        ax = twoVec_hist_compare(seq(i,osc_ids),seq(i,n_osc_ids),...
            seqOrder{1,i},100,...
            {'osc','n-osc'},'plot');
    end
end

end