function plot_fit_over_gen(obj,whichFit,gen_num)
% this function plots the maximum fitness in every generation

% 'whichFit' - row vector of the fitnesess to plot.
% 'gen_num' - up to which generation number to plot

x_data = 1:gen_num;
y_data = zeros(4,gen_num);

for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    
    for j=1:numel(obj.data_names)
        % get the fit for all generations:
        y = obj.data{1,j}.GA.Fit(:,FitNum,x_data);
        genes = obj.data{1,j}.GA.Seqs(:,:,x_data);
        
        % get wether this fit should be minimize/maximize:
            % '1'-fot max;      '-1'- for min
        minORmax = obj.data{1,j}.GA.FitMinMax;
        
        % get y-axis data:
        switch minORmax(1,j)
            case 1
                [bestFit,bestID] = max(squeeze(y),[],1);
                for g=1:length(x_data)
                    bestGenes(g,:) = genes(bestID(g),:,g);
                end
                y_data(j,:) = bestFit;
            case -1
                [M,I] = min(y,[],1);
                y_data(j,:) = squeeze();
        end
        
%         % OPTIONAL: check best gene passing
%         checkGenerationContinuaty(obj,x_data,genes,bestGenes);
    end
    
    figure; hold on;
    for j=1:numel(obj.data_names)
        plot(x_data, y_data(j,:));
    end
    grid minor;
    legend(obj.Legends, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Fitness');
    title(['Fitness "',obj.fitnessOrder{1,FitNum},'" over Generation']);
    set(gca,'FontSize',12, 'FontWeight','bold');

end

    function checkGenerationContinuaty(obj,x_data,genes,bestGenes)
        % check that the best gene is passed to the next geeration
        % check if the best gene is there in the next generation:
        for g=2:length(x_data)
            bsxfun(@eq,squeeze(genes(:,:,g)),bestGenes(g-1,:),bestFit);
        end
        
        % plot each graph next to a table with values
        figure;
        T = array2table([x_data',bestGenes,bestFit'],'VariableNames',...
            {'genNum','tau','b','c','NR','a','kc','kt','bestFit'});
        uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
            'RowName',T.Properties.RowNames,'Units', 'Normalized',...
            'Position',[0, 0, 1, 1]);
        
        figure;
        plot(x_data, bestFit); grid minor;
        xlabel('Generation');   ylabel('Fitness');
        title(['Fitness "',obj.fitnessOrder{1,FitNum},'" over Generation']);
        set(gca,'FontSize',12, 'FontWeight','bold');
    end

end

