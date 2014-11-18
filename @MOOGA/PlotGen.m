function PlotGen( GA, varargin ) %#ok<INUSL>
%PLOTGEN Prepares plots for a specific genome
%   Prepares plots for the genome based on the analysis performed.
%   Possible plots: 'IC','LC','EigV','EVLocus'

% Plot formatting
AxesFont = 12;
LabelFont = 16;
TitleFont = 20;
LineWidth = 3;
LineStyles = {'-','--',':','-.','-*','-o'};
Markers = {'+','o','d','^','v'};
Colors = {[0 0 1],[1 0 0],[0 0.7 0.8],[0 0.7 0],[0.8 0 0.8]};
Legends = {'\theta_1','\theta_2','d\theta_1/dt',...
                    'd\theta_2/dt','\phi_C_P_G'};

Data = varargin{1};
if nargin<3
    PlotAll();
else
    NPlots = nargin/2-1;
    if mod(NPlots,1) ~= 0
        error('Wrong number of arguments provided');
    end
    
    for arg = 1:NPlots
        PlotByName(varargin{2*arg},varargin{2*arg+1});
    end
end

    function PlotAll()
        PlotIC();
        PlotLC();
    end

    function PlotByName(Name,Options)
        switch lower(Name)
            case {'ic','initcond','initial conditions'}
                PlotIC(Options);
            case {'lc','limitcycle','limit cycle'}
                PlotLC(Options);
        end
    end

    function PlotIC(Options)
        if nargin<1
            % Default options
            sub = 0;
        else
            op = 1;
            while op<=length(Options)
                switch lower(Options{op})
                    case {'sub','subplot'}
                        % Read next option for the subplot values
                        sub = 1;
                        subIDs = Options{op+1};
                        op = op+1;
                        NsubP = length(unique(subIDs));
                end
                op = op+1;
            end
        end
        
        figure
        hold on
        
        IC = Data.IC;
        Slopes = Data.Slopes;
        
        Ncoords = size(IC,1)/max(Data.Period(:,1));
        h = zeros(Ncoords,1);
        for p = 1:length(Data.Zones)
            pZone = Data.Zones{p};
            for z = 1:length(pZone)
                Coords = pZone{z};
                for c = 1:Ncoords
                    StCoord = (p-1)*Ncoords+c;
                    if sub == 1
                        % Plot each coordinate in a subplot
                        subplot(NsubP,1,subIDs(c));
                        hold on
                    end
                    h(c) = plot(Slopes(Coords),IC(StCoord,Coords),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
            end
        end
        
        % Stretch the plot and add legends + formatting
        if sub == 1
            for sp = 1:NsubP
                subplot(NsubP,1,sp);
                IDs = find(subIDs == sp);
                axis([min(Slopes) max(Slopes) ylim])
                legend(h(IDs),Legends(IDs),'FontSize',AxesFont);
                
                set(gca,'FontSize',AxesFont);
                if sp<NsubP
                    set(gca,'XTickLabel',[]);
                else
                    xlabel('Slope [deg]','FontSize',LabelFont);
                end
            end
        else
            axis([min(Slopes) max(Slopes) ylim])
            legend(h,Legends,'FontSize',AxesFont);
            set(gca,'FontSize',AxesFont);
            xlabel('Slope [deg]','FontSize',LabelFont);
        end
    end

    function PlotLC(Options)
    end
        
end



