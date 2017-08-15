function plot_Tend_ratio(data,gen_num,titleAdd)
% plot the ratio of: 
%               T(end)          the time in which the sim stopped
%   ratio = --------------- = -------------------------------------
%               Sim.Tend        the run sim max time (as defined)
% 
% Inputs:
% *) 'data' - cell array contain MOGa results
% *) 'gen_num' - the generation to focus
% *) 'titleAdd' - the title for each case
% 
% 
% 

% generation numbers:
X_data = 1:gen_num;

%% plot the maxRatio + mean and stdev
figure;
for i=1:4
    Y = data{1,i}.GA.Tend_ratio(:,1,X_data);
    
    y_data_max = squeeze(max(Y,[],1));
    y_data_mean = squeeze(mean(Y,1));
    y_data_std = squeeze(std(Y,[],1));
    
    ax = subplot(2,2,i); hold on;
    plot(ax,X_data,y_data_max);
    plot(ax,X_data,y_data_mean);
    errorbar(ax,X_data,y_data_mean,y_data_std);
    legend('Max','mean');
    title(['T(end)/Sim.Tend  for case: ',titleAdd{1,i}]);
    xlabel('generation num');
    ylabel('T(end)/Sim.Tend');
    grid minor
    
end

% waht counts as a successfull Sim: (CB walked)
good_ratio_threshold = 0.99;


%% ALSO: plot the percentage of sims that the ratio is close to '1':
figure; hold on;
ax = gca;
for i=1:4
    Y = squeeze(data{1,i}.GA.Tend_ratio(:,1,X_data));
    
    good_ids = (Y > good_ratio_threshold);
    amount_of_good_ids = sum(good_ids,1);
    percent_good = 100*(amount_of_good_ids/size(Y,1));
    
    plot(ax,X_data,percent_good);
  
end
title(['percent of Sims that got: T(end)/Sim.Tend > ',...
    num2str(good_ratio_threshold),' in each generation']);
xlabel('generation num');
ylabel('%');
grid minor
legend(titleAdd,'Location','Best');
hold off;

%% LAST: plot the ratio distribution in each generation (on heatmap):
figure; 

for i=1:4
    Y = squeeze(data{1,i}.GA.Tend_ratio(:,1,X_data));
    
    Prob = distribution_over_genNum(Y,1);
    
    subplot(2,2,i); hold on
    imagesc(Prob');
    title({'the distribution of the ratio: T(end)/Sim.Tend over generation',...
        ['for case: ',titleAdd{1,i}]});
    xlabel('generation num');
    ylabel('T(end)/Sim.Tend');
    grid minor
    axis([0,gen_num,0,100]);
    hold off
end

%% LAST 2: plot the ratio distribution in each generation(bar graph):
figure; 

for i=1:4
    Y = squeeze(data{1,i}.GA.Tend_ratio(:,1,X_data));
    
    Prob = distribution_over_genNum(Y,1);

    subplot(2,2,i); hold on
    bar3(Prob');
    view(-56,43);
    title({'the distribution oaf the ratio: T(end)/Sim.Tend over generation',...
        ['for case: ',titleAdd{1,i}]});
    xlabel('generation num');
    ylabel('T(end)/Sim.Tend [%]');
    zlabel('probability');
    grid minor
    axis([0,gen_num,0,100]);
    hold off
end


end

