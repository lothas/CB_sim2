function my_MoE_plot_test_perf(obj,expertCount,NNoutput,NNtargets,cluster_i_train_ind,gateOut,competetiveFlag)
% plot only graphs which are relevant to test sessions:

switch competetiveFlag
    case {1,2}
        % plot the regression graph and color each expert's cluster in a different color 
        figure; hold on
        for j=1:expertCount
            if length(cluster_i_train_ind{1,j}) > 0
                out_temp = NNoutput(:,cluster_i_train_ind{1,j});
                targ_temp = NNtargets(:,cluster_i_train_ind{1,j});

                h = plot(targ_temp,out_temp,'Color',obj.colors(j,:),'LineStyle','none');
                h.Marker = 'o';
            end
        end
        hold off;
        xlabel('targets'); ylabel('ouput'); legend(obj.legendNames);
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
        legend(obj.legendNames);

        % figure with different filled color for each "dominant expert
        % (based on "g")
        [g_max,g_max_ind] = max(gateOut,[],1);
        figure; hold on
        for j=1:expertCount
            for i=1:size(NNoutput,2)
                if g_max_ind(1,i) == j
                    if g_max(1,i) > 0.5
                        plot(NNtargets(1,i),NNoutput(1,i),'k-o','MarkerFaceColor',obj.colors(j,:));
                    else
                        plot(NNtargets(1,i),NNoutput(1,i),'k-o');
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

