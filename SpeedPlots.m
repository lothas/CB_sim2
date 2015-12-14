function SpeedPlots()
%     data = load('GA_11_24_10_25.mat');
%     ids = 1:19;
%     ids(2) = [];
    
    data = load('VGA_12_14_12_02.mat');
    ids = 1:18;
    GA = data.GA;
    
    % Get first Pareto front
    Data = GA.Fit(:,1:2,GA.Progress);
    Data = [Data (1:size(Data,1))'];
    Fronts = GA.Pareto(Data);
    
    Best = [Fronts{1};Fronts{2}];
    N = length(Best);
    Fits = GA.Fit(Best,1:2,GA.Progress);
    
    Top = [Fits,Best];
    sTop = sortrows(Top,1);
    
    % Take out points that are too close
    speed_min = min(Fits(:,1));
    speed_max = max(Fits(:,1));
    exp_delta = 0.5*(speed_max - speed_min)/N;
    
    i = 1;
    while i<size(sTop,1)
        dist = sqrt((sTop(i,1)-sTop(i+1,1))^2+(sTop(i,2)-sTop(i+1,2))^2);
        if dist < exp_delta
            sTop(i+1,:) = [];
        else
            i = i+1;
        end
    end
    Seqs = GA.Seqs(sTop(:,end),ids,GA.Progress);
    
    % Plot parameters (without feedback)
    figure
    splots = {1:3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    ylabels = {'f_osc','T_a','p_a','dp_a','T_h1','p_h1','dp_h1','T_h2','p_h2','dp_h2'};
    x = sTop(:,1);
    for i = 1:10
        subplot(4,3,splots{i});
        y = Seqs(:,i);
        yLPF = LPF(y);
        scatter(x,y,'.');
        hold on
        plot(x(1:length(yLPF)),yLPF,'LineWidth',2);
        ylabel(ylabels{i});
    end
        
    function filt = LPF(sig)
        NF = 5;
        filt = sig(1:length(sig)-NF);
        for s = 1:length(sig)-NF
            filt(s) = sum(sig(s:s+NF-1))/NF;
        end
    end
%     scatter(Top(:,1),Top(:,2))
%     hold on
%     scatter(sTop(:,1),sTop(:,2),'*')
        
%     scatter(Fits(:,1),Fits(:,2))
end