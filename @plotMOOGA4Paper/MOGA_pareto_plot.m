function Fronts = MOGA_pareto_plot(obj,fit1Num,fit2Num,gen_num,plotType)
% plot pareto fronts

% Inputs:
% *) 'data' - GA results in cell array
% *) 'fit1Num','fit2Num' - the fitness Id numbers
% *) 'fitnessOrder' - the order in which the fitnesses are organized
% *) 'gen_num' - the relevant generation number
% *) 'plotType' : 
%               #) 'showAll' - show specific fronts for every case in 
%                               a different subplot.
%               #) 'showBest' - show specific front for every case on the 
%                               same figure (each color for different case)
% 
% Outputs:
% *) the Pareto Fronts. Cell array containf the genes Ids

LnWidth = 4;
FontSizeT = 14;
FontSize = 12;
FontSize2 = 10;

figure();

switch plotType
    case {'showAll','subplot'}
        for i = 1:4

            subplot(2,2,i);
            GA = obj.data{1,i}.GA;

            % % % % uncomment this if you want to normalized the fitnesses by maxFit.
            % % Normalize fitness
            % GA.Fit(:,fit1Num,:) = ...
            %     GA.Fit(:,fit1Num,:)/max(max(GA.Fit(:,fit1Num,:)));
            % GA.Fit(:,fit2Num,:) = ...
            %     GA.Fit(:,fit2Num,:)/max(max(GA.Fit(:,fit2Num,:)));

            % % % % % dont forget to add the correct label later:)
            % xlabel(['normalzed ',fitnessOrder{1,fit1Num}],'FontSize',FontSize)
            % ylabel(['normalzed ',fitnessOrder{1,fit2Num}],'FontSize',FontSize)

            Fi = [1, 10, 20];
            NF = length(Fi);

            % labels generation:
            for f = 1:NF
                Label{1,f} = ['front #',num2str(Fi(f))];
            end


            % figure('units','normalized','Position',[0.0953, 0.3870, 0.3161, 0.3139]);
            Data = [GA.Fit(:,[fit1Num,fit2Num],gen_num),(1:GA.Population)'];
            Fronts = GA.Pareto(Data);

            hold on
            for f = 1:NF
                FrData = sortrows(Data(Fronts{Fi(f)},:));
                x = FrData(:,1);
                y = FrData(:,2);
                Color = f/NF*[1, 0, 0] + (1-f/NF)*[0, 0.8, 0];
                plot(x,y,'-x','Color',Color,'LineWidth',3,'MarkerSize',12);
            end
            % % Add the selected genome
            % SelectedGen = 4;
            % GenFit = GA.Fit(SelectedGen,[fit1Num,fit2Num],gen_num);
            % plot(GenFit(1),GenFit(2),'ko','MarkerSize',12,...
            %     'MarkerFaceColor',[0 0 0]);
            xlabel(obj.fitnessOrder{1,fit1Num},'FontSize',FontSize)
            ylabel(obj.fitnessOrder{1,fit2Num},'FontSize',FontSize)
            set(gca,'FontSize',FontSize2,'LineWidth',LnWidth/2)
            title({['Pareto plot of generation #',num2str(gen_num)],...
                ['for case #',obj.titleAdd{1,i}]},...
                'FontSize',FontSizeT);
            legend(Label,'Location','Best');

        end
    case {'showBest'}
        Color = {[0,0.447,0.741],...
            [0.850,0.325,0.098],...
            [0.929,0.694,0.125],...
            [0.494,0.184,0.556]};
        Legend = {'GA only','GA + NN','GA + rescale','GA + NN + rescale'};
        for i = 1:4
            
            GA = obj.data{1,i}.GA;

            %Fi = [1,2,3];
            Fi = [4];
            NF = length(Fi);

            Data = [GA.Fit(:,[fit1Num,fit2Num],gen_num),(1:GA.Population)'];
            Fronts = GA.Pareto(Data);

            hold on
            x = [];
            y = [];
            for f = 1:NF
                FrData = sortrows(Data(Fronts{Fi(f)},:));
                x = [x;FrData(:,1)];
                y = [y;FrData(:,2)];
            end
            plot(x,y,'x','Color',Color{1,i},'MarkerSize',12);
            
        end
        xlabel(obj.fitnessOrder{1,fit1Num},'FontSize',FontSize)
        ylabel(obj.fitnessOrder{1,fit2Num},'FontSize',FontSize)
        set(gca,'FontSize',FontSize2,'LineWidth',LnWidth/2)
        title({['Pareto plot of generation #',num2str(gen_num)],...
            ['Front num #',num2str(Fi)]},...
            'FontSize',FontSizeT);
        legend(Legend,'Location','Best');
        
    otherwise
        error('invalid plotType');
end



end

