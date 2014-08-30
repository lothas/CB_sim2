function [ qData ] = Quantile(GA,Data,DesQ)
%QUANTILE Returns the high-valued part of the Data representing the DesQ quantile
%   For example: If DesQ is 0.2, Quantile will return the top 20% of Data
%   NOTE: Data should be provided a unique ID as the last column.
Method = 3;
% 1 - Leave top quantile by eliminating low values from each category
% 2 - Leave top quantile by eliminating inverse Pareto fronts

if nargin<3
    DesQ = 0.5;
end

if nargin<2
    close all
    
    NRand = 1000;
    Nr = 5;

    % Generate random 3D data
    Data = zeros(NRand,3);
    Ns = NRand/Nr;
    for r=1:Nr
        theta = pi/2*rand(Ns,1);
%         phi = pi/2-pi*rand(Ns,1);
        phi = pi/2*rand(Ns,1);
%         theta = pi/2*((0.5+rand(Ns,1)).^2-0.25)/2;
%         phi = pi/2*((0.5+rand(Ns,1)).^2-0.25)/2;
        Data(1+Ns*(r-1):Ns*r,1) = r*cos(phi).*cos(theta);
        Data(1+Ns*(r-1):Ns*r,2) = r*cos(phi).*sin(theta);
        Data(1+Ns*(r-1):Ns*r,3) = r*sin(phi);
    end
%     Data = max(0,min(1./rand(NRand,3),10)-5);
    x = max(0,1/5-rand(NRand,1))*5;
    y = max(0,1/5-rand(NRand,1))*5;
    z = max(0,rand(NRand,1)+1./(1+x+y)-1.5);
    Data = [x y z];

    % Give unique ID to each sample
    Data = [Data (1:size(Data,1))'];

    figure
    subplot(1,3,1);
    scatter3(Data(:,1),Data(:,2),Data(:,3),'k');
    view(130,54)
    hold on
end

% Find the quantile that selects only the best individuals
% from each category
switch Method
    case 1
        if isempty(DesQ)
            % Worst case scenario will eliminate 30% of the data
            % (if each column is "independent")
            quant = 0.3/size(Data,2);
        else
            options = optimset('MaxIter',10,'TolFun',1e-3);
            quant0 = DesQ;
            quant = fzero(@quantf,quant0,options,Data,DesQ);
        end

        [qData] = GetQuant(Data,quant);
    case 2
        % Get Inverse fronts
        Fronts = GA.Pareto(Data,1);

        % Remove fronts until enough IDs have been taken
        IDs = zeros(size(Data,1));
        f = 1;
        nIDs = 0;
        IDsReq = floor((1-DesQ)*size(Data,1));
        while nIDs<IDsReq
            FL = length(Fronts{f});
            IDs(1+nIDs:nIDs+FL) = Fronts{f};
            nIDs = nIDs+FL;
            f = f+1;
        end
        IDs = IDs(1:IDsReq);
        qData = Data;
        qData(IDs,:) = [];
    case 3
        DesOut = DesQ*size(Data,1);
            
        % Check number of values above zero
        zData = Data>0;
        
        % Count zeros
        zData_count = sum(zData(:,1:end-1),2); 
        
        Nc = size(zData,2)-1; % Number of columns (fits)
        Nabove = zeros(Nc+1,1);
        IDabove = cell(Nc+1,1);
        
        NOut = 0;
        IDsOut = [];
        % Number of data points to be deleted:
        NOutRequired = floor((1-DesQ)*size(Data,1));
        groupsd = 0; % groups deleted
        
        for co = 1:Nc+1
            % How many scored above 0 "co-1" times?
            IDabove{co} = find(zData_count==co-1);
            Nabove(co) = length(IDabove{co});
            NOut = NOut+length(IDabove{co});
            if NOut<NOutRequired
                IDsOut = [IDsOut;IDabove{co}]; %#ok<AGROW>
                groupsd = groupsd+1;
            end
        end
        
%         if groupsd == 0
%             % IDs with all zero values are already too many
%             % We'll duplicate the good result IDs
%             GoodIDs = cell2mat(IDabove(2:end));
%             rep = floor(DesOut/length(GoodIDs));
%             IDs = repmat(GoodIDs,rep,1);
%             IDs = [IDs; randsample(GoodIDs,DesOut-length(IDs))];
%             qData = Data(IDs,:);
%             return
%         end
            
        % Delete IDs with too many zero values, but save the "best" ones
        qData = Data;
        if ~isempty(IDsOut)
            Front = GA.Pareto(Data(IDsOut,:),0,1); % Get one front
            IDsOut = setdiff(IDsOut,Front{1});
        end
        
        if groupsd<Nc+1
            % The next group has too many IDs in it, so we'll pick out the
            % best and only remove the worst
            Fronts = GA.Pareto(Data(IDabove{groupsd+1},:));
            Nkeep = (size(Data,1)-length(IDsOut))-DesOut;
            
            f = 1;
            IDsOut2 = IDabove{groupsd+1};
            while 1
                FL = length(Fronts{f});
                if length(IDsOut2)-FL<Nkeep
                    break
                end
                IDsOut2 = setdiff(IDsOut2, Fronts{f});
                f = f+1;
            end
            
            % Select at random from last front
            LastIDs = randsample(Fronts{f},length(IDsOut2)-Nkeep);
            IDsOut = [IDsOut; setdiff(IDsOut2, LastIDs)];
        end
        
        qData(IDsOut,:) = [];
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
        Qdiff = ResQ-DesQ;
        if Qdiff<0
            Qdiff = 1;
        end
    end
        
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
%         tri = delaunay(x,y);
%         trisurf(tri, x, y, z,'FaceColor',[1-i/N; i/2/N; i/N]);
        scatter3(x,y,z,'*','MarkerEdgeColor',[1; i/N; i/N]);
        view(130,54)
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
%         tri = delaunay(x,y);
%         trisurf(tri, x, y, z,'FaceColor',[1-i/N; i/2/N; i/N]);
        scatter3(x,y,z,'*','MarkerEdgeColor',[1; i/N; i/N]);
        view(130,54)
    end
    
    mean(qData(:,1:3))-mean(Data(:,1:3))    
end

end

