function my_MoE_plot_test_perf(obj,expertCount,NNoutput,NNtargets,cluster_i_train_ind,gateOut,competetiveFlag)
% plot only graphs which are relevant to training sessions:

switch expertCount
    case {2,3} % in case of small number of expert, make colors clear:
        colors = [1,0,0;0,1,0;0,0,1];
    otherwise
        colors = rand(expertCount,3);
end

legendNames = cell(1,expertCount);
for j=1:expertCount
    legendNames{1,j} = ['#',num2str(j),' expert'];
end

switch competetiveFlag
    case {1,2}
        % plot the regression graph and color each expert's cluster in a different color 
        figure; hold on
        for j=1:expertCount
            out_temp = NNoutput(:,cluster_i_train_ind{1,j});
            targ_temp = NNtargets(:,cluster_i_train_ind{1,j});

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

        % figure with different filled color for each "dominant expert
        % (based on "g")
        [g_max,g_max_ind] = max(gateOut,[],1);
        colors = rand(expertCount,3);
        figure; hold on
        for j=1:expertCount
            for i=1:size(NNoutput,2)
                if g_max_ind(1,i) == j
                    if g_max(1,i) > 0.5
                        plot(NNoutput(1,i),NNtargets(1,i),'k-o','MarkerFaceColor',colors(j,:));
                    else
                        plot(NNoutput(1,i),NNtargets(1,i),'k-o');
                    end
                end
            end
        end
        xlabel('target');    ylabel('output'); grid on;
        title({'regression graph with different color for every dominant expert';'empty circle mean g<0.5'});
        hold off

    otherwise
        error('wrong "competetiveFlag", try again');   
end


end

