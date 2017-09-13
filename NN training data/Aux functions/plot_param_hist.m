function plot_param_hist(seq,periods,seqOrder)
% PLOT_PARAM_HIST plot the histograms of all of the parameters

N = length(seqOrder);

figure;
for i=1:(N+1)
    subplot(4,4,i);
    if (i > N)
        histogram(periods,100,'Normalization','pdf');
        title(['histogram of: ','periods']);
    else
        histogram(seq(i,:),100,'Normalization','pdf');
        title(['histogram of: ',seqOrder{1,i}]);
    end
end

end

