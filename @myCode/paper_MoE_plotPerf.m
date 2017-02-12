function paper_MoE_plotPerf(obj,graphKeys)
% plot the train and/or test performance of the papers function

% fetcching data from class:
sampl_test = obj.sampl_test;
targ_test = obj.targ_test;
targ_train = obj.targ_train;

numOfIteretions = obj.numOfIteretions;
expertCount = obj.expertCount;
Moe_perf_over_iter = obj.paper_MoE_out.Moe_perf_over_iter;
netOut_train = obj.paper_MoE_out.out_from_train;
netOut_test = obj.paper_MoE_out.out_from_test;
gateOut_train = obj.paper_MoE_out.gateOut_from_train;
gateOut_test = obj.paper_MoE_out.gateOut_from_test;

% set legend and x-axis for plot:
iterNum = 1:numOfIteretions;
legendNames = cell(1,expertCount);
for j=1:expertCount
    legendNames{1,j} = ['#',num2str(j),' expert'];
end

% plot the requested graphs:
for i=1:length(graphKeys)
    switch graphKeys{1,i}
        case {'MSE_over_iter','MSE_train'}
            % MSE over iteration num:
            figure;
            plot(iterNum,Moe_perf_over_iter,'b-o');
            xlabel('#iteration');   ylabel('MSE err');
            title('MSE error over #iteration');
            
        case {'regGraph','reg'}
            % regression graph with gateOut bar graph
            figure;
            subplot(2,2,1);
            h1 = plot(targ_test,netOut_test,'LineStyle','none'); grid minor;
            h1.Marker = 'o';
            xlabel('targets'); ylabel('outputs');
            title('regression graph: targets over output for Testing');

            subplot(2,2,3);
            bar(gateOut_test','stacked'); xlabel('#sample');
            legend(legendNames);
            title('gate output over sample num - test');

            subplot(2,2,2);
            h2 = plot(targ_train,netOut_train,'LineStyle','none'); grid minor;
            h2.Marker = 'o';
            xlabel('targets'); ylabel('outputs');
            title('regression graph: targets over output for Training');

            subplot(2,2,4);
            bar(gateOut_train','stacked'); xlabel('#sample');
            legend(legendNames);
            title('gate output over sample num - train');
            
        case {'reg_graph_from_NNtoolbox'}
            % regression as plotted from NN toolbox
                figure;
                plotregression(targ_train,netOut_train,'train',...
                    targ_test,netOut_test,'test');
                
        case {'one_dimentional_problem'}
            % if we have one dimentional problem (such as the 'V' problem)
            figure; 
            h2 = plot(sampl_test,targ_test,'LineStyle','none'); hold on;
            h2.Marker = 'o';
            h3 = plot(sampl_test,netOut_test,'Color','r','LineStyle','none');
            h3.Marker = 'x';
            legend('targets','output');
            
        otherwise
            disp('the options are:');
            disp('"MSE_over_iter", and/or');
            disp('"regGraph", and/or');
            disp('"reg_graph_from_NNtoolbox", and/or');
            disp('"one_dimentional_problem", (only if number of input = one)');
            error('illegal input key!');
    end

end
end

