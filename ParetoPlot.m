clear all
close all

LnWidth = 3;
MrkrSize = 12;
FontSizeT = 24;
FontSize = 20;
FontSize2 = 16;

load('GA_11_13_19_31_St1');
Fi = [1, 10, 20];
NF = length(Fi);

figure
Data = [GA.Fit(:,[1 2],GA.Progress),(1:GA.Population)'];
Fronts = GA.Pareto(Data);
subplot(1,2,1);
hold on
for f = 1:NF
    FrData = sortrows(Data(Fronts{Fi(f)},:));
    x = FrData(:,1);
    y = FrData(:,2);
    Color = f/NF*[1, 0, 0] + (1-f/NF)*[0, 0.8, 0];
    plot(x,y,'-x','Color',Color,'LineWidth',3,'MarkerSize',12);
end
xlabel('Velocity fitness','FontSize',FontSize)
ylabel('Energy efficiency fitness','FontSize',FontSize)
set(gca,'FontSize',FontSize2)
title('A','FontSize',FontSizeT);

Fi = [1, 4, 15];
Data = [GA.Fit(:,[4 5],GA.Progress),(1:GA.Population)'];
Fronts = GA.Pareto(Data);
subplot(1,2,2);
hold on
for f = 1:NF
    FrData = sortrows(Data(Fronts{Fi(f)},:));
    x = FrData(:,1);
    y = FrData(:,2);
    Color = f/NF*[1, 0, 0] + (1-f/NF)*[0, 0.8, 0];
    plot(x,y,'-x','Color',Color,'LineWidth',3,'MarkerSize',12);
end
xlabel('Uphill fitness','FontSize',FontSize)
ylabel('Downhill fitness','FontSize',FontSize)
set(gca,'FontSize',FontSize2)
title('B','FontSize',FontSizeT);

% Data = [GA.Fit(:,[7 8],GA.Progress),(1:GA.Population)'];
% Fronts = GA.Pareto(Data);
% subplot(1,3,3);
% hold on
% for f = 1:NF
%     FrData = sortrows(Data(Fronts{Fi(f)},:));
%     x = FrData(:,1);
%     y = FrData(:,2);
%     Color = f/NF*[1, 0, 0] + (1-f/NF)*[0, 0.8, 0];
%     plot(x,y,'-x','Color',Color,'LineWidth',3,'MarkerSize',12);
% end
% xlabel('Downhill slope fitness','FontSize',FontSize)
% ylabel('ZMP fitness','FontSize',FontSize)
% set(gca,'FontSize',FontSize2)
% title('C','FontSize',FontSizeT);