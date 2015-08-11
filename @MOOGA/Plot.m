function Plot( GA,varargin )
%PLOT Outputs different plots of the algorithm's performance

% Plots format
AxesFont = 16;
LabelFont = 18;
TitleFont = 20; %#ok<NASGU>
LineWidth = 4;
LineStyles = {'-','--',':','-.','-*','-o'};
Markers = {'+','o','d','^','v'};
Colors = {[0 0 1],[1 0 0],[0 0.7 0.8],[0 0.7 0],[0.8 0 0.8]};
Legends = {'\theta_1','\theta_2','d\theta_1/dt',...
                    'd\theta_2/dt','\phi_C_P_G'};
mylegends = {'F_{Vel}','F_{Energy}','F_{Conv}','F_{Slope}'};

if nargin<2
    PlotFit(0);
else
    switch varargin{1}
        case {'Fit','fit','fitness'}
            PlotFit(1);
        case {'FitMax','fitmax'}
            PlotFit(0);
        case {'Gen','gen'}
            if nargin<3
                error('No specific genome provided');
            else
                if length(varargin{2}) == 1
                    ID = varargin{2}(1);
                    Gen = GA.Progress;
                else
                    ID = varargin{2}(1);
                    Gen = varargin{2}(2);
                end

                % Open analysis data if it exists
                Data = GA.Analyze(Gen,ID);

                if nargin == 3
                    GA.PlotGen(Data);
                else
                    GA.PlotGen(Data,varargin{3:end});
                end
            end
    end
end

    function PlotFit(Type)
        % Type 0: Just max
        % Type 1: Max and mean
        figure('units','normalized','Position',[0.1, 0.1, 0.45, 0.6])
        hold on
        Generations = 1:GA.Generations;
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

                if fi == 1
                    Thish = plot(Generations,FitMaxi,'Color',Colors{f},...
                        'LineWidth',LineWidth);
                end
                h(f) = Thish(1);
                
                if Type == 1
                    if fi>1 
                        plot(Generations,FitMaxi,'-.','Color',Colors{f},...
                            'LineWidth',LineWidth);
                    end
                    plot(Generations,FitMeani,'--','Color',Colors{f},...
                        'LineWidth',LineWidth);
                end

                if Max>1 && (fi == 1 || Type == 1)
                    % Add a text balloon with the max value
                    ht = text(Generations(maxID)-0.2,0.98,num2str(Max,'%.2f'),...
                        'HorizontalAlignment','center',...
                        'BackgroundColor',[1 1 1],'FontSize',12,...
                        'VerticalAlignment','top');
                    set(ht,'FontSize',12);
                end
            end
            
%             legends{f} = MOOGA.GetFitFcnName(GA.FitFcn{f,2});
        end
        legend(h,mylegends,'Location','SouthEast','FontSize',LabelFont);
        xlabel('Generations','FontSize',LabelFont);
        ylabel('Normalized Fitness','FontSize',LabelFont);
        set(gca,'FontSize',AxesFont,'LineWidth',2);
    end
        
end

