function Plot( GA,which )
%PLOT Outputs different plots of the algorithm's performance
switch which
    case {'Fit','fit','fitness'}
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
            
            MaxMax = max(max(FitMax));
            if MaxMax>1
                % Normalize values larger than 1
                maxID = find(FitMax==MaxMax,1,'first');
                row = 1+mod(maxID,size(FitMax,1));
                col = ceil(maxID/size(FitMax,1));
                FitMax(row,:) = FitMax(row,:)/MaxMax;
            end
                
            Thish = plot(Generations,FitMax,'Color',Colors{f},...
                'LineWidth',2);
            h(f) = Thish(1);
            plot(Generations,FitMean,'--','Color',Colors{f},...
                'LineWidth',2);
            
            if MaxMax>1
                % Add a text balloon with the max value
                text(Generations(col)-0.2,0.98,num2str(MaxMax,'%.2f'),...
                    'HorizontalAlignment','center',...
                    'BackgroundColor',[1 1 1],'FontSize',12,...
                    'VerticalAlignment','top');
            end
            
            legends{f} = MOOGA.GetFitFcnName(GA.FitFcn{f,2});
        end
        legend(h,legends,'Location','NorthWest');
end

end

