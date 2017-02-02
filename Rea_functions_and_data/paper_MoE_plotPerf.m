function paper_MoE_plotPerf(NNinput,NNtargets,NNoutput,gateOut,errsOverIter)

numOfIteretions = length(errsOverIter);

expertCount = size(gateOut,1);
legendNames = cell(1,expertCount);
for j=1:expertCount
    legendNames{1,j} = ['#',num2str(j),' expert'];
end

figure;
plot(1:numOfIteretions,errsOverIter,'b-o');
xlabel('#iteration');   ylabel('MSE err');
title('MSE error over #iteration');

figure;
subplot(2,1,1);
h = plot(NNtargets,NNoutput,'LineStyle','none'); grid minor;
h.Marker = 'o';
xlabel('targets'); ylabel('outputs');
title('regression graph: Targets over NNoutputs');

subplot(2,1,2);
bar(gateOut','stacked'); xlabel('#sample');
legend(legendNames);

[~,~] = NN_perf_calc(NNtargets,NNoutput',1,1);

if size(NNinput,1) == 1
    figure; 
    h2 = plot(NNinput,NNtargets,'LineStyle','none'); hold on;
    h2.Marker = 'o';
    h3 = plot(NNinput,NNoutput,'Color','r','LineStyle','none');
    h3.Marker = 'x';
    legend('targets','output');
end

end

