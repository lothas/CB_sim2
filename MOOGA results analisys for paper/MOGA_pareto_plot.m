function Fronts = MOGA_pareto_plot(data,fit1Num,fit2Num,fitnessOrder,...
    gen_num,titleAdd)
% plot pareto fronts

% Inputs:
% *) 'data' - GA results in cell array
% *) 'fit1Num','fit2Num' - the fitness Id numbers
% *) 'fitnessOrder' - the order in which the fitnesses are organized
% *) 'gen_num' - the relevant generation number
% *)'titleAdd' - for titles cases names

% Outputs:
% *) the Pareto Fronts. Cell array containf the genes Ids

LnWidth = 4;
FontSizeT = 14;
FontSize = 12;
FontSize2 = 10;

figure();

for i = 1:4
    
    subplot(2,2,i);
    GA = data{1,i}.GA;

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
    xlabel(fitnessOrder{1,fit1Num},'FontSize',FontSize)
    ylabel(fitnessOrder{1,fit2Num},'FontSize',FontSize)
    set(gca,'FontSize',FontSize2,'LineWidth',LnWidth/2)
    title({['Pareto plot of generation #',num2str(gen_num)],...
        ['for case #',titleAdd{1,i}]},...
        'FontSize',FontSizeT);
    legend(Label,'Location','Best');
    
end


end

