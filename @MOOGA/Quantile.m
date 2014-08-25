function [ qData ] = Quantile(GA,Data,DesQ)
%QUANTILE Returns the high-valued part of the Data representing the DesQ quantile
%   For example: If DesQ is 0.2, Quantile will return the top 20% of Data
%   NOTE: Data should be provided a unique ID as the last column.

if nargin<3
    DesQ = 0.1;
end

if nargin<2
    close all
    
    NRand = 5000;
    Nr = 5;

    % Generate random 3D data
    Data = zeros(NRand,3);
    Ns = NRand/Nr;
    for r=1:Nr
        theta = pi/2*rand(Ns,1);
        phi = pi/2*rand(Ns,1);
%         theta = pi/2*((0.5+rand(Ns,1)).^2-0.25)/2;
%         phi = pi/2*((0.5+rand(Ns,1)).^2-0.25)/2;
        Data(1+Ns*(r-1):Ns*r,1) = r*cos(phi).*cos(theta);
        Data(1+Ns*(r-1):Ns*r,2) = r*cos(phi).*sin(theta);
        Data(1+Ns*(r-1):Ns*r,3) = r*sin(phi);
    end
%     Data = min(1./rand(NRand,3),10);

    % Give unique ID to each sample
    Data = [Data (1:size(Data,1))'];

    figure
    subplot(1,3,1);
    scatter3(Data(:,1),Data(:,2),Data(:,3),'k');
    hold on
end

% Find the quantile that selects only the best individuals
% from each category
if isempty(DesQ)
    % Worst case scenario will eliminate 30% of the data
    % (if each column is "independent")
    quant = 0.3/size(Data,2);
else
    options = optimset('MaxIter',10,'TolFun',1e-3);
    quant0 = DesQ;
    quant = fzero(@quantf,quant0,options,Data,DesQ);
end

    function [qData] = GetQuant(Data,quant)
        qIDs = 1:size(Data,1);
        for c = 1:size(Data,2)-1
            qnt = quantile(Data(:,c),quant);
            qIDs = intersect(qIDs,find(Data(:,c)>=qnt));
        end
        qData = Data(qIDs,:);
    end

    function [Qdiff] = quantf(quant,Data,DesQ)
        quant = min(max(quant,0),1);
        [qData] = GetQuant(Data,quant);
        if ~isempty(qData)
            ResQ = size(qData,1)/size(Data,1);
        else
            ResQ = 0;
        end
        Qdiff = DesQ-ResQ;
    end
        
[qData] = GetQuant(Data,quant);

if nargin<2
    scatter3(qData(:,1),qData(:,2),qData(:,3),'g*');
    
    % Give unique ID to each sample
    Data = [Data (1:size(Data,1))'];
    
    Fronts = GA.Pareto(Data);
    N=length(Fronts);

    subplot(1,3,2);
    hold on
    for i=1:N
        if length(Fronts{i})<3
            continue;
        end
        x = Data(Fronts{i},1);
        y = Data(Fronts{i},2);
        z = Data(Fronts{i},3);
        tri = delaunay(x,y);
        trisurf(tri, x, y, z,'FaceColor',[1-i/N; i/2/N; i/N]);
    end
    
    % Give unique ID to each sample
    qData = [qData (1:size(qData,1))'];
    
    Fronts = GA.Pareto(qData);
    N=length(Fronts);

    subplot(1,3,3);
    hold on
    for i=1:N
        if length(Fronts{i})<3
            continue;
        end
        x = qData(Fronts{i},1);
        y = qData(Fronts{i},2);
        z = qData(Fronts{i},3);
        tri = delaunay(x,y);
        trisurf(tri, x, y, z,'FaceColor',[1-i/N; i/2/N; i/N]);
    end
end

end

