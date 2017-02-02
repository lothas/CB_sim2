function my_MoE_plotPerf(NNoutput,NNtargets,gateOut,...
    cluster_i_train_ind,Moe_perf_over_iter,gateNN_perf_vec,testOrTrain,competetiveFlag )

% inputs:
% 1) testOrTrain - plot graphs from training or test sessions

expertCount = size(gateOut,1);

colors = rand(expertCount,3);
legendNames = cell(1,expertCount);
for j=1:expertCount
    legendNames{1,j} = ['#',num2str(j),' expert'];
end

% plot only graphs which are relevant to training sessions:
if strcmp(testOrTrain,'test') || strcmp(testOrTrain,'both')
    switch competetiveFlag
        case {1,2}
            % plot the regression graph and color each expert's cluster in a different color 
            figure; hold on
            for j=1:expertCount
                lengthGroup_j = size(cluster_i_train_ind{1,j},2);
                out_temp = NNoutput(:,1:lengthGroup_j);
                targ_temp = NNtargets(:,1:lengthGroup_j);

                h = plot(targ_temp,out_temp,'Color',colors(j,:),'LineStyle','none');
                h.Marker = 'o';
            end
            hold off;
            xlabel('targets'); ylabel('ouput'); legend(legendNames);
            title('regression graph: Targets over NNoutputs');

        case 3
            figure;
            subplot(2,1,1);
            h = plot(NNtargets,NNoutput,'LineStyle','none'); grid minor;
            h.Marker = 'o';
            xlabel('targets'); ylabel('outputs');
            title('regression graph: Targets over NNoutputs');

            subplot(2,1,2);
            bar(gateOut','stacked'); xlabel('#sample');
            legend(legendNames);

        otherwise
            error('wrong "competetiveFlag", try again');   
    end
end

if strcmp(testOrTrain,'train') || strcmp(testOrTrain,'both')
    
    numOfIteretions = length(Moe_perf_over_iter); 
    switch competetiveFlag
        case {1,2}
            figure;
            subplot(2,2,1)
            plot(1:numOfIteretions,expert_i_GroupSize,'-o'); hold on;
            title('cluster size over #iteration');
            xlabel('#iteretion');   ylabel('group size [#points]');
            legend(expertsNames);

            subplot(2,2,2) % expert perf over #interation
            plot(1:numOfIteretions,Experts_perf_mat,'-o'); hold on;
            title('expert MSE over #interation');
            xlabel('#iteretion');   ylabel('performance [MSE]');

            subplot(2,2,3) % expert perf over #interation
            plot(1:numOfIteretions,double(emptyGroupIndecator),'-o'); hold on;
            title('indication on empty clusters: "1" means empty');
            xlabel('#iteretion');   ylabel('"1"=empty   "0"-not empty');

    end
    
    figure; 
    subplot(2,1,1) % total MSE error over #iteration
    plot(1:numOfIteretions,Moe_perf_over_iter,'b-o');
    xlabel('#iteration'); ylabel('MoE MSE error');
    title('total MSE error over #iteration');
    
    subplot(2,1,2) % gateNet perf over #interation
    plot(1:numOfIteretions,gateNN_perf_vec,'-o'); hold on;
    title('gateNet perf (crossEntropy) over #interation');
    xlabel('#iteretion');   ylabel('performance [crossentropy]');
end

[~,~] = NN_perf_calc(NNtargets,NNoutput,1,1);

end

