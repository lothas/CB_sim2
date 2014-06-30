function Plot( GA,which )
%PLOT Outputs different plots of the algorithm's performance
switch which
    case {'Fit','fit','fitness'}
        figure()
        hold on
        Generations = 1:GA.Generations;
        Colors = {[1 0 0],[0 0 1],[0 0 0],[0 0.7 0],[0.6 0.4 0]};
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

    function splitstr = strsplit(str,delimiters)
        Nd = length(delimiters);
        places = 0;
        for d = 1:Nd
            pos = find(str == delimiters{d});
            places = [places, pos]; %#ok<AGROW>
        end
        places = sort(places);
        
        if places(end)<length(str)
            splitstr = cell(1,length(places));
            for p = 1:length(places)-1
                splitstr{p} = str(places(p)+1:places(p+1)-1);
            end
            splitstr{p+1} = str(places(end)+1:end);
        else
            splitstr = cell(1,length(places)+1);
            for p = 1:length(places)-1
                splitstr{p} = str(places(p)+1:places(p+1)-1);
            end
        end
    end

end

