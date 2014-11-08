function Plot( GA,which )
%PLOT Outputs different plots of the algorithm's performance
switch which
    case {'Fit','fit','fitness'}
        PlotFit(1);
    case {'FitMax','fitmax'}
        PlotFit(0);
end

    function PlotFit(Type)
        % Type 0: Just max
        % Type 1: Max and mean
        figure()
        hold on
        Generations = 1:GA.Generations;
        Colors = {[1 0 0],[0 0 1],[0 0 0],[0 0.7 0],...
                [0.6 0.4 0],[0 0.4 0.6],[0 0.7 0.3]};
        h = zeros(1,GA.NFit);
        legends = cell(1,GA.NFit);
        for f = 1:GA.NFit
            FitInd = GA.FitFcn{f,1};
            LF = length(FitInd);
            FitMax = reshape(max(GA.Fit(:,FitInd,:)),LF,[]);
            FitMean = reshape(mean(GA.Fit(:,FitInd,:)),LF,[]);
            
            for fi = 1:LF
                FitMaxi = FitMax(fi,:);
                FitMeani = FitMean(fi,:);
                Max = max(FitMaxi);
                if Max>1
                    % Normalize values larger than 1
                    maxID = find(FitMaxi==Max,1,'first');
                    FitMaxi = FitMaxi/Max;
                    FitMeani = FitMeani/Max;
                end

                Thish = plot(Generations,FitMaxi,'Color',Colors{f},...
                    'LineWidth',2);
                h(f) = Thish(1);
                if Type == 1
                    plot(Generations,FitMeani,'--','Color',Colors{f},...
                        'LineWidth',2);
                end

                if Max>1
                    % Add a text balloon with the max value
                    text(Generations(maxID)-0.2,0.98,num2str(Max,'%.2f'),...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[1 1 1],'FontSize',12,...
                        'VerticalAlignment','top');
                end
            end
            
            legends{f} = MOOGA.GetFitFcnName(GA.FitFcn{f,2});
        end
        legend(h,legends,'Location','NorthWest');
    end
end

