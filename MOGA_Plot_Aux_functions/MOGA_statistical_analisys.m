function [ output_args ] = ...
    MOGA_statistical_analisys(seqOrder,MML,...
    resultsFileNames,categoryLegend,whichGraph)
%MOGA_STATISTICAL_ANALISYS get names of MOGA run files devider into
%categories and perform statistical analisys on the results
% 
% INPUTS:
% *) 'MML' - matsouka analisys class
% *) 'seqOrder' - 1Xn cell array with the parameters order in the sequence
% *) 'resultsFileNames' - 1Xm cell array contain the names of the results
%                       files.
%           for example:
%               resultsFileNames = {GA_res,NN_res};
%               GA_res = {'GAfile1.mat','GAfile2.mat',...};
%               NN_res = {'NNfile1.mat','NNfile2.mat',...};
% *) 'whichGraph' - which graphs to show
%      #)'mean_of_max_fitness' = extract the max fitness of each generation
%           from every GA file and then take the mean for each generation.
%      #)
% 
% *) 'categoryLegend' - 1Xm cell array contain the names of categories
% 
% 
% 
% 
% 

% Set colors:
COLORS = [242,74,38 ; 53,110,251; 223,174,21] ./255;
COLORS_light = [249,195,174 ; 191,208,251; 249,225,149] ./255;

% set default options
set(0,'defaultlinelinewidth',2);

% % define the class for CPG simulation:
% MML = MatsuokaML();
% MML.perLim = [0.68 0.78];
% MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
% MML.tStep = 0.05;
% MML.tEnd = 15;
% seqOrder = {'tau' ,'b', 'c', 'w1', 'w2', 'w3', 'w4','k_tau','k_{c}'};

% % Load the files:
categoryNum = numel(resultsFileNames);
MOGA_results = cell(1,categoryNum);
for i=1:categoryNum
    Legends = {'GA'};
    MOGA_results{1,i} = ...
        plotMOOGA4Paper(MML,resultsFileNames{1,i},Legends,seqOrder);
end

% % get the results:


% get x-axis data:
last_gen = 20; % last generation to check

% fitNum = 1:3 = {'velFit','NrgFit','rangeVelFit'}
Title = {'velFit','NrgFit','rangeVelFit'};
whichFit2Plot = 1:3;

% generations num vector:
genNum = 1:last_gen;

for j=1:length(whichFit2Plot)
    meanFit = zeros(categoryNum,last_gen);
    stdFit = zeros(categoryNum,last_gen);
    switch whichGraph
        % graphs of Fitness comparison:
        case 'mean_of_max_fitness_errBars' % % Plot mean of maximum fitness for all MOGA runs
                for i=1:categoryNum
                    [~,maxFit] = ...
                        MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);
                    % calc mean and stdev:
                    meanFit(i,:) = mean(maxFit);
                    stdFit(i,:) = std(maxFit);
                    maxFitVec{1,i} = maxFit;
                end

                figure; hold on;
                y_max_lim = plot_data(categoryNum,maxFitVec,meanFit,stdFit,'errBars');

                xlabel('Generation number');
                ylabel('Fitness avg');
                legend(categoryLegend);
                grid minor
                axis([0.5,last_gen,0,y_max_lim]);
                set(gca,'fontsize',13);
                title(['\fontsize{12}',Title{1,j},': max fitness avarage across multiple runs']);
                hold off;

        case 'mean_of_max_fitness_thinLines' % Plot mean of maximum fitness for all MOGA runs:
                    for i=1:categoryNum
                        [~,maxFit] = ...
                            MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);

                        % calc mean and stdev:
                        meanFit(i,:) = mean(maxFit);
                        stdFit(i,:) = std(maxFit);
                        maxFitVec{1,i} = maxFit;
                    end

                    figure; hold on;
                    y_max_lim = plot_data(categoryNum,maxFitVec,meanFit,stdFit,'thinLines');

                    xlabel('Generation number');
                    ylabel('Fitness avg');
                    % legend([h{1,1},h{1,2}],categoryLegend);
                    legend(categoryLegend);
                    grid minor
                    axis([0.5,last_gen,0,y_max_lim]);
                    set(gca,'fontsize',13);
                    title({['\fontsize{12}',Title{1,j},': max fitness avarage across multiple runs'],...
                        '\fontsize{11} Bold lines = avg      clear lines = each sim'});
                    hold off;

        case 'mean_of_mean_fitness_errBars' % mean of meanFit at each generation
            % % Plot mean of "TopPop Mean Fitness" (TPMF) fitness for all MOGA runs:
                for i=1:categoryNum
                    [~,AvgFit] = ...
                        MOGA_results{1,i}.plot_mean_fit_over_gen(whichFit2Plot(j),last_gen,'top_pop');

                    % calc mean and stdev:
                    meanFit(i,:) = mean(AvgFit);
                    stdFit(i,:) = std(AvgFit);
                    avgFitVec{1,i} = AvgFit;
                end

                figure; hold on;
                y_max_lim = plot_data(categoryNum,avgFitVec,meanFit,stdFit,'errBars');

                xlabel('Generation number');
                ylabel('Fitness Mean');
                legend(categoryLegend);
                grid minor
                axis([1,last_gen,0,y_max_lim]);
                set(gca,'fontsize',13);
                title(['\fontsize{12}',Title{1,j},'_{mean} of the avarage fitness']);
                hold off;

        case 'mean_of_max_fitness_Signrank_test_per_Generation' % add statistical test
                for i=1:categoryNum
                    [~,maxFit] = ...
                        MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);

                    % calc mean and stdev:
                    meanFit(i,:) = mean(maxFit);
                    stdFit(i,:) = std(maxFit);

                    maxFitVec{1,i} = maxFit;
                end

                % perf a SignRank test between the first two:
                [prob,rejectFlag] = wilcox_sign_rank_test(maxFitVec,last_gen);

                figure; hold on;
                y_max_lim = plot_data(categoryNum,maxFitVec,meanFit,stdFit,'errBars');
                
                % place Markers when the results are statistically
                % significant:
                pointHight = y_max_lim+0.1;
                place_signifi_markers(last_gen,pointHight,rejectFlag);
                
                xlabel('Generation number', 'FontSize', 13);
                ylabel('Fitness Mean', 'FontSize', 13);
                legend(categoryLegend);
                set(gca,'fontsize',13);
                title({['\fontsize{13}',Title{1,j},': max fitness avarage across multiple runs'],...
                    ['\fontsize{10} green = statistically significant ;     red = not stat. sig']},'interpreter','tex');
                grid minor
                axis([0.5,last_gen,0,pointHight+0.05]);
                
                % % Uncomment if you want to print the p-values:
%                 place_pValue_text(last_gen,pointHight,prob);
                hold off;

        case 'mean_of_max_fitness_Signrank_test_All_Generation' % perform the signRankTest to all generation combained
                for i=1:categoryNum
                    [~,maxFit] = ...
                        MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);

                    % calc mean and stdev:
                    meanFit(i,:) = mean(maxFit);
                    stdFit(i,:) = std(maxFit);

                    maxFitVec{1,i} = maxFit;
                end

                % perf a SignRank test between the first two:
                disp('Note: performing test only on the 1st two GA runs');
                V1 = [];
                V2 = [];
                for k=1:last_gen
                    V1 = [V1;maxFitVec{1,1}(:,k)];
                    V2 = [V2;maxFitVec{1,2}(:,k)];
                end
                [p,h] = signrank(V1,V2);

                disp(['at fit: ',Title{1,j}]);
                if h == 1
                    disp('    The null hypothesis was regected');
                    disp('    the results are statistically significant!');
                else
                    disp('    The null hypothesis was NOT regected');
                    disp('    the results are NOT statistically significant!');
                end
    end
end

mean_runTime = zeros(1,last_gen);
std_runTime = zeros(1,last_gen);
switch whichGraph % graphs about the runTime of the MOGA
    case 'mean_of_max_runTime_Signrank_test_per_Generation' % % graphs of RunTime comparison
            for i=1:categoryNum
                [~,runTime] = ...
                    MOGA_results{1,i}.plot_gen_runTime(last_gen);
                % calc mean and stdev:
                mean_runTime(i,:) = mean(runTime);
                std_runTime(i,:) = std(runTime);
                
                runTimeVec{1,i} = runTime;
            end

            % perf a SignRank test between the first two:
            [prob,rejectFlag] = wilcox_sign_rank_test(runTimeVec,last_gen);
    
            figure; hold on;
            y_max_lim = plot_data(categoryNum,runTimeVec,...
                mean_runTime,std_runTime,'errBars')

            % place Markers when the results are statistically
            % significant:
            pointHight = y_max_lim*1.000001;
            place_signifi_markers(last_gen,pointHight,rejectFlag);

            xlabel('Generation number', 'FontSize', 13);
            ylabel('RunTime [sec]', 'FontSize', 13);
            legend(categoryLegend);
            set(gca,'fontsize',13);
            title('\fontsize{13} avg. runTime of each generation');
            grid minor
            axis([0.5,last_gen,0,y_max_lim]);

            % % Uncomment if you want to print the p-values:
%                 place_pValue_text(last_gen,pointHight,prob);
            hold off;
        
    case 'mean_of_max_runTime_thinLines' % % graphs of RunTime comparison
            for i=1:categoryNum
                [~,runTime] = ...
                    MOGA_results{1,i}.plot_gen_runTime(last_gen);
                % calc mean and stdev:
                mean_runTime(i,:) = mean(runTime);
                std_runTime(i,:) = std(runTime);
                
                runTimeVec{1,i} = runTime;
            end

            figure; hold on;
            y_max_lim = plot_data(categoryNum,runTimeVec,...
                mean_runTime,std_runTime,'thinLines')

            xlabel('Generation number', 'FontSize', 13);
            ylabel('RunTime [sec]', 'FontSize', 13);
            legend(categoryLegend);
            set(gca,'fontsize',13);
            title('\fontsize{13} avg. runTime of each generation');
            grid minor
            axis([0.5,last_gen,0,y_max_lim]);
            hold off;
            
    case 'mean_runTime_until_the_ith_gen'
        % % calc the run time of the GA until the i-th generation
        % % i.e. it includes the time of the previuos generations
        
        I = ones(last_gen,last_gen);
        I = triu(I); % create an upper triangular matrix
        
        for i=1:categoryNum
            [~,runTime] = ...
                MOGA_results{1,i}.plot_gen_runTime(last_gen);
            
            runTime_tot = runTime*I;
            
            % calc mean and stdev:
            mean_runTime(i,:) = mean(runTime_tot);
            std_runTime(i,:) = std(runTime_tot);

            runTimeVec{1,i} = runTime_tot;
        end

        % perf a SignRank test between the first two:
        [prob,rejectFlag] = wilcox_sign_rank_test(runTimeVec,last_gen);

        figure; hold on;
        y_max_lim = plot_data(categoryNum,runTimeVec,...
            mean_runTime,std_runTime,'errBars')

        % place Markers when the results are statistically
        % significant:
        pointHight = y_max_lim*1.000001;
        place_signifi_markers(last_gen,pointHight,rejectFlag);

        xlabel('Generation number', 'FontSize', 13);
        ylabel('RunTime [sec]', 'FontSize', 13);
        legend(categoryLegend);
        set(gca,'fontsize',13);
        title('\fontsize{13} avg. runTime until the i^{th} generation');
        grid minor
        axis([0.5,last_gen,0,y_max_lim]);

        % % Uncomment if you want to print the p-values:
%                 place_pValue_text(last_gen,pointHight,prob);
        hold off;
end

    function [pValue,rejectFlag] = wilcox_sign_rank_test(res_vec,last_gen)
        % perf a SignRank test between the first two:
        
        % % [p,h] = signrank(___) also returns a logical value indicating 
        % %     the test decision. 
        % %
        % % 'p': p-value of a paired, two-sided test for the null hypothesis
        % %     that x – y comes from a distribution with zero median.
        % %
        % %  h = 1 indicates a rejection of the null hypothesis, and
        % %  h = 0 indicates a failure to reject the null hypothesis at the
        % %      5% (percent) significance level.
        % %  You can use any of the input arguments in the previous syntaxes.

        disp('Note: performing test only on the 1st two GA runs');
        for g=1:last_gen
            [pValue(1,g),rejectFlag(1,g)] =...
                signrank(res_vec{1,1}(:,g),res_vec{1,2}(:,g));
        end
    end

    function [h] = place_signifi_markers(last_gen,pointHight,rejectFlag)
        for gen=1:last_gen
                    if rejectFlag(1,gen)==1
                        h = scatter(gen,pointHight,'filled','Marker','o');
                        h.MarkerFaceColor = [124,252,80]./255;
                    else
                        h = scatter(gen,pointHight,'filled','Marker','o');
                        h.MarkerFaceColor = [230,51,7]./255;
                    end
                    h.MarkerEdgeColor = zeros(1,3);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
        end
    end

    function place_pValue_text(last_gen,pointHight,pValue)
        text(0.55,pointHight + 0.04,'p-values:','FontSize',12);
        
        textHight = 0.02;
        for gen=1:last_gen
            % place text of the p-value
            text(gen,...
                pointHight - textHight,...
                sprintf('%.3f',pValue(1,gen)),...
                'FontSize',11,...
                'HorizontalAlignment','center');
            textHight = textHight*(-1);
        end
    end

    function y_max_lim = plot_data(categoryNum,res_vec,Mean,stdev,plotType)
        for c=1:categoryNum
            plot(genNum,Mean(c,:),'Color',COLORS(c,:));
            switch plotType
                case 'errBars' % plot with error bars
                    h3 = errorbar(genNum,Mean(c,:),stdev(c,:));
                    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                case 'thinLines' % plot the max fits in thin lines
                    for fit=1:size(res_vec{1,c},1)
                        h3 = plot(genNum,res_vec{1,c}(fit,:),...
                            'Color',COLORS_light(c,:),'LineWidth',2);
                        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    end
            end
        end
        
        % find the maximum point:
        y_limits = get(gca,'YLim');
        y_max_lim = y_limits(1,2);
    end


if false
switch whichGraph
    % graphs of Fitness comparison:
    case 'mean_of_max_fitness_errBars'
        Title = {'velFit','NrgFit','rangeVelFit'};
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:length(whichFit2Plot)

            meanFit = zeros(categoryNum,last_gen);
            stdFit = zeros(categoryNum,last_gen);
            
            for i=1:categoryNum
                [~,maxFit] = ...
                    MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);
                
                % calc mean and stdev:
                meanFit(i,:) = mean(maxFit);
                stdFit(i,:) = std(maxFit);
            end

            figure; hold on;
            for i=1:categoryNum
                plot(genNum,meanFit(i,:),'Color',COLORS(i,:));
                h1 = errorbar(genNum,meanFit(i,:),stdFit(i,:));
                set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end

            xlabel('Generation number');
            ylabel('Fitness Mean');
            legend(categoryLegend);
            title([Title{1,j},'_{mean}']);
            grid minor
            axis([1,last_gen,0,0.4]);
            set(gca,'fontsize',13);
            hold off;
        end
    case 'mean_of_max_fitness_thinLines'
            Title = {'velFit','NrgFit','rangeVelFit'};
            genNum = 1:last_gen;
            % % Plot mean of maximum fitness for all MOGA runs:
            for j=1:length(whichFit2Plot)   
                meanFit = zeros(categoryNum,last_gen);
                stdFit = zeros(categoryNum,last_gen);

                for i=1:categoryNum
                    [~,maxFit] = ...
                        MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);

                    % calc mean and stdev:
                    meanFit(i,:) = mean(maxFit);
                    stdFit(i,:) = std(maxFit);
                    
                    maxFitcell{1,i} = maxFit;
                end
                
                figure; hold on;
                for i=1:categoryNum
                    plot(genNum,meanFit(i,:),'Color',COLORS(i,:));
                    % errorbar(genNum,meanFit(i,:),stdFit(i,:));
                    
                    % plot the max fits in thin lines:
                    for k=1:size(maxFitcell{1,i},1)
                        h1 = plot(genNum,maxFitcell{1,i}(k,:),...
                            'Color',COLORS_light(i,:),'LineWidth',2);
                        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    end
                end

                xlabel('Generation number');
                ylabel('Fitness Mean');
                % legend([h{1,1},h{1,2}],categoryLegend);
                legend(categoryLegend);
                title([Title{1,j},'_{mean}']);
                grid minor
                axis([1,last_gen,0,0.4]);
                set(gca,'fontsize',13);
                hold off;
            end
    case 'mean_of_mean_fitness_errBars' % mean of meanFit at each generation
        % % Plot mean of "TopPop Mean Fitness" (TPMF) fitness for all MOGA runs:
        Title = {'velFit','NrgFit','rangeVelFit'};
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:length(whichFit2Plot)
            meanFit = zeros(categoryNum,last_gen);
            stdFit = zeros(categoryNum,last_gen);
            
            for i=1:categoryNum
                [~,AvgFit] = ...
                    MOGA_results{1,i}.plot_mean_fit_over_gen(whichFit2Plot(j),last_gen,'top_pop');
                
                % calc mean and stdev:
                meanFit(i,:) = mean(AvgFit);
                stdFit(i,:) = std(AvgFit);
            end

            figure; hold on;
            for i=1:categoryNum
                plot(genNum,meanFit(i,:),'Color',COLORS(i,:));
                h1 = errorbar(genNum,meanFit(i,:),stdFit(i,:));
                set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end

            xlabel('Generation number');
            ylabel('Fitness Mean');
            legend(categoryLegend);
            title([Title{1,j},'_{mean} of the avarage fitness']);
            grid minor
            axis([1,last_gen,0,0.4]);
            set(gca,'fontsize',13);
            hold off;
        end
    case 'mean_of_max_fitness_Signrank_test_per_Generation'
        Title = {'velFit','NrgFit','rangeVelFit'};
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:length(whichFit2Plot)
            meanFit = zeros(categoryNum,last_gen);
            stdFit = zeros(categoryNum,last_gen);
            
            for i=1:categoryNum
                [~,maxFit] = ...
                    MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);
                
                % calc mean and stdev:
                meanFit(i,:) = mean(maxFit);
                stdFit(i,:) = std(maxFit);
                
                maxFitVec{1,i} = maxFit;
            end
            
            % perf a SignRank test between the first two:
            disp('Note: performing test only on the 1st two GA runs');
            for k=1:last_gen
                [p(1,k),h(1,k)] =...
                    signrank(maxFitVec{1,1}(:,k),maxFitVec{1,2}(:,k));
            end
            
            figure; hold on;
            for i=1:categoryNum
                
                plot(genNum,meanFit(i,:),'Color',COLORS(i,:));
                h1 = errorbar(genNum,meanFit(i,:),stdFit(i,:));
                set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                % place Markers when the results are statistically
                % significant:
                for k=1:last_gen
                    if h(1,k)==1
                        h2 = scatter(k,0.5,'filled','Marker','o');
                        h2.MarkerFaceColor = [124,252,80]./255;
                    else
                        h2 = scatter(k,0.5,'filled','Marker','o');
                        h2.MarkerFaceColor = [230,51,7]./255;
                    end
                    h2.MarkerEdgeColor = zeros(1,3);
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
            end

            xlabel('Generation number');
            ylabel('Fitness Mean');
            legend(categoryLegend);
            title([Title{1,j},'_{mean}']);
            grid minor
            axis([1,last_gen,0,0.51]);
            set(gca,'fontsize',13);
            hold off;
        end    
    case 'mean_of_max_fitness_Signrank_test_All_Generation'
        Title = {'velFit','NrgFit','rangeVelFit'};
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:length(whichFit2Plot)
            meanFit = zeros(categoryNum,last_gen);
            stdFit = zeros(categoryNum,last_gen);
            
            for i=1:categoryNum
                [~,maxFit] = ...
                    MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot(j),last_gen);
                
                % calc mean and stdev:
                meanFit(i,:) = mean(maxFit);
                stdFit(i,:) = std(maxFit);
                
                maxFitVec{1,i} = maxFit;
            end
            
            % perf a SignRank test between the first two:
            disp('Note: performing test only on the 1st two GA runs');
            V1 = [];
            V2 = [];
            for k=1:last_gen
                V1 = [V1;maxFitVec{1,1}(:,k)];
                V2 = [V2;maxFitVec{1,2}(:,k)];
            end
            [p,h] = signrank(V1,V2);
            
            disp(['at fit: ',Title{1,j}]);
            if h == 1
                disp('    The null hypothesis was regected');
                disp('    the results are statistically significant!');
            else
                disp('    The null hypothesis was NOT regected');
                disp('    the results are NOT statistically significant!');
            end
        end        
end
end

end

