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
            FitMax = reshape(max(GA.Fit(:,f,:)),1,[]);
            FitMean = reshape(mean(GA.Fit(:,f,:)),1,[]);
            h(f) = plot(Generations,FitMax,'Color',Colors{f},...
                'LineWidth',2);
            plot(Generations,FitMean,'--','Color',Colors{f},...
                'LineWidth',2);
            FitFcn = strsplit(func2str(GA.FitFcn{f}),{'.','('});
            legends{f} = FitFcn{3};
        end
        legend(h,legends,'Location','NorthWest');
end

end

