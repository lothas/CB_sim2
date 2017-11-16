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
last_gen = 10;


switch whichGraph
    case 'mean_of_max_fitness_errBars'
        Title = {'velFit','NrgFit','rangeVelFit'};
        last_gen = 10; % last generation to check
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:3   % fitNum = 1:3 = {'velFit','NrgFit','rangeVelFit'}
            whichFit2Plot = j;
            
            meanFit = zeros(categoryNum,10);
            stdFit = zeros(categoryNum,10);
            
            for i=1:categoryNum
                [~,maxFit] = ...
                    MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot,last_gen);
                
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
            axis([1,10,0,0.4]);
            set(gca,'fontsize',13);
            hold off;
        end
    case 'mean_of_max_fitness_thinLines'
            Title = {'velFit','NrgFit','rangeVelFit'};
            last_gen = 10; % last generation to check
            genNum = 1:last_gen;
            % % Plot mean of maximum fitness for all MOGA runs:
            for j=1:3   % fitNum = 1:3 = {'velFit','NrgFit','rangeVelFit'}
                whichFit2Plot = j;

                meanFit = zeros(categoryNum,10);
                stdFit = zeros(categoryNum,10);

                for i=1:categoryNum
                    [~,maxFit] = ...
                        MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot,last_gen);

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
                axis([1,10,0,0.4]);
                set(gca,'fontsize',13);
                hold off;
            end
    case 'mean_of_mean_fitness_errBars' % mean of meanFit at each generation
        % % Plot mean of "TopPop Mean Fitness" (TPMF) fitness for all MOGA runs:
        Title = {'velFit','NrgFit','rangeVelFit'};
        last_gen = 10; % last generation to check
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:3   % fitNum = 1:3 = {'velFit','NrgFit','rangeVelFit'}
            whichFit2Plot = j;
            
            meanFit = zeros(categoryNum,10);
            stdFit = zeros(categoryNum,10);
            
            for i=1:categoryNum
                [~,AvgFit] = ...
                    MOGA_results{1,i}.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'top_pop');
                
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
            axis([1,10,0,0.4]);
            set(gca,'fontsize',13);
            hold off;
        end
    case 'mean_of_max_fitness_Signrank_test_per_Generation'
        Title = {'velFit','NrgFit','rangeVelFit'};
        last_gen = 10; % last generation to check
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:3   % fitNum = 1:3 = {'velFit','NrgFit','rangeVelFit'}
            whichFit2Plot = j;
            
            meanFit = zeros(categoryNum,10);
            stdFit = zeros(categoryNum,10);
            
            for i=1:categoryNum
                [~,maxFit] = ...
                    MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot,last_gen);
                
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
            axis([1,10,0,0.51]);
            set(gca,'fontsize',13);
            hold off;
        end    
    case 'mean_of_max_fitness_Signrank_test_All_Generation'
        Title = {'velFit','NrgFit','rangeVelFit'};
        last_gen = 10; % last generation to check
        genNum = 1:last_gen;
        % % Plot mean of maximum fitness for all MOGA runs:
        for j=1:3   % fitNum = 1:3 = {'velFit','NrgFit','rangeVelFit'}
            whichFit2Plot = j;
            
            meanFit = zeros(categoryNum,10);
            stdFit = zeros(categoryNum,10);
            
            for i=1:categoryNum
                [~,maxFit] = ...
                    MOGA_results{1,i}.plot_fit_over_gen(whichFit2Plot,last_gen);
                
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



% if false % % Plot the AVG runTime of each generation
%     [x_data,GA_only_runTime] = GA_only.plot_gen_runTime(last_gen);
%     [~,GA_NN_runTime] = GA_NN.plot_gen_runTime(last_gen);
% 
%     GA_only_mean = mean(GA_only_runTime);
%     GA_only_std = std(GA_only_runTime);
% 
%     GA_NN_mean = mean(GA_NN_runTime);
%     GA_NN_std = std(GA_NN_runTime);
% 
%     figure;
%     plot(x_data,GA_only_mean); hold on;
%     errorbar(x_data,GA_only_mean,GA_only_std);
%     plot(x_data,GA_NN_mean);
%     errorbar(x_data,GA_NN_mean,GA_NN_std);
%     xlabel('Generation number');
%     ylabel('Avg runTime [sec]');
%     legend('MOGA','MOGA with NN assist');
%     title('Avg runTime');
%     grid minor
%     axis([1,10,0,inf]);
%     set(gca,'fontsize',13)
% 
%     clear x_data GA_only_runTime GA_NN_runTime 
%     clear GA_only_mean GA_only_std GA_NN_mean GA_NN_std
% end
% 
% if false % % Plot the runTime of until the i-th generation
%     I = ones(last_gen,last_gen);
%     I = triu(I); % create an upper triangular matrix
%     
%     [x_data,GA_only_runTime] = GA_only.plot_gen_runTime(last_gen);
%     [~,GA_NN_runTime] = GA_NN.plot_gen_runTime(last_gen);
%     
%     % summing each genTime with the previuos time to get the time until
%     % that point (Time Elapesd)
%     GA_only_elapsed_time = GA_only_runTime*I;
%     GA_only_mean = mean(GA_only_elapsed_time);
%     GA_only_std = std(GA_only_elapsed_time);
%     
%     GA_NN_elapsed_time = GA_NN_runTime*I;
%     GA_NN_mean = mean(GA_NN_elapsed_time);
%     GA_NN_std = std(GA_NN_elapsed_time);
% 
%     figure;
%     plot(x_data,GA_only_mean); hold on;
%     errorbar(x_data,GA_only_mean,GA_only_std);
%     plot(x_data,GA_NN_mean);
%     errorbar(x_data,GA_NN_mean,GA_NN_std);
%     xlabel('Generation number');
%     ylabel('Avg runTime [sec]');
%     legend('MOGA','MOGA with NN assist');
%     title('Avg runTime');
%     grid minor
%     axis([1,10,0,inf]);
%     set(gca,'fontsize',13)
% 
%     clear x_data GA_only_runTime GA_NN_runTime 
%     clear GA_only_elapsed_time GA_NN_elapsed_time I
%     clear GA_only_mean GA_only_std GA_NN_mean GA_NN_std
% end

end

