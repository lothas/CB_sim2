function [ TopIDs, Weights ] = GetTopPop( GA, N )
%GETTOPPOP Receives a population's fitness and returns the N top IDs
%   GetTopPop divides the given genomes into Pareto fronts based
%   on their fitness. It then selects the top N genomes from the
%   layers one by one. If the layer is larger than required,
%   genomes that excell at each fitness individually are picked.

    if GA.NFit == 1
        FitID = [GA.Fit(:,GA.FitIDs,GA.Progress), (1:GA.Population)'];
        Sorted = sortrows(FitID,-1);
        TopIDs = Sorted(1:N,2);
        Weights = ones(N,1);
    else
        TopIDs = zeros(GA.Population,1);
        Weights = ones(GA.Population,1);
        
        Data = GA.Fit(:,GA.FitIDs,GA.Progress);
        % Give unique ID to each sample
        Data = [Data (1:size(Data,1))'];
        
        f = 1;
        IDindex = 1;
        switch GA.JOAT % Jack of all trades
            %%%%%%%%%%%%%%%%% Specialization "allowed" %%%%%%%%%%%%%%%%%
            case 0
                % Divide population into pareto fronts
                Fronts = GA.Pareto(Data);
                while IDindex<=N
                    % Take IDs from pareto fronts until TopPop is reached
                    [f,TopIDs,Weights,IDindex] = ...
                        GetFront(Fronts,f,TopIDs,Weights,IDindex);
                end

                if f==2 && IDindex>N
                    % The first pareto front had more than N individuals
                    % Select the best genomes from each fitness
                    TopIDs = GetByBest(GA,TopIDs,IDindex,N);
                end
            %%%%%% Specialization allowed for single "elite" stratum %%%%%%
            case 1
                % Divide population into pareto fronts
                Fronts = GA.Pareto(Data);
                
                % Take the first front
                [~,TopIDs,Weights,IDindex] = ...
                    GetFront(Fronts,1,TopIDs,Weights,IDindex);
                Weights = Weights - 1;
                
                Lf1 = length(Fronts{1});
                if Lf1<=N
                    % Remove those IDs from Data
                    Data(Fronts{1},:) = [];
                    
                    % Get top quantile
                    Jacks = GA.Quantile(Data,GA.Quant);
                    
                    % Divide population into pareto fronts
                    Fronts = GA.Pareto(Jacks);
                    while IDindex<=N
                        % Take IDs from pareto fronts until TopPop is reached
                        [f,TopIDs,Weights,IDindex] = ...
                            GetFront(Fronts,f,TopIDs,Weights,IDindex);
                    end
                    Weights = Weights + 1;                
                else
                    % First front will fill the quota
                    TopIDs = GetByBest(GA,TopIDs,IDindex,N);
                end
            %%%%%%%%%%%%%% Jacks of all trades are enforced %%%%%%%%%%%%%%
            case 2
                % Get top quantile
                Jacks = GA.Quantile(Data,GA.Quant);

                % Divide population into pareto fronts
                Fronts = GA.Pareto(Jacks);
                while IDindex<=N
                    % Take IDs from pareto fronts until TopPop is reached
                    [f,TopIDs,Weights,IDindex] = ...
                        GetFront(Fronts,f,TopIDs,Weights,IDindex);
                end
                
                if f==2 && IDindex>N
                    % The first pareto front had more than N individuals
                    % Select the best genomes from each fitness
                    TopIDs = GetByBest(GA,TopIDs,IDindex,N);
                end
        end
        
        TopIDs = TopIDs(1:N);
        % Invert the weights (first fronts get higher weight)
        Weights = f-Weights;
        Weights = Weights(1:N);
    end
    
    %%%%%%%%%%%%%%%%%% GETFRONT %%%%%%%%%%%%%%%%%%
    function [f,TopIDs,Weights,IDindex] = ...
        GetFront(Fronts,f,TopIDs,Weights,IDindex)
        % Take IDs from pareto front
        res = length(Fronts{f});
        TopIDs(IDindex:IDindex+res-1) = Fronts{f}(randperm(res));
        Weights(IDindex:IDindex+res-1) = repmat(f,res,1);
        % Fronts are randomly permutated to remove biases
        IDindex = IDindex+res;
        f = f+1;
    end
    
    %%%%%%%%%%%%%%%%%% GETBYBEST %%%%%%%%%%%%%%%%%%
    function TopTopIDs = GetByBest(GA,TopIDs,IDindex,N)
        IDindex2 = 1;
        FitIndex = 1;
        TopTopIDs = zeros(N,1);
        
    %         if min(Fit(1,TopIDs(1:IDindex-1)))>0.95
    %             % Don't select genomes by height fitness (Fit(1))
    %             Fit0 = 2;
    %         else
                Fit0 = 1;
    %         end
    
        while IDindex2<=N
            % Select the highest genomes for each fitness
            FitInd = GA.FitIDs(FitIndex);
            
            ThisID = find(GA.Fit(:,FitInd,GA.Progress) == ...
                max(GA.Fit(TopIDs(1:IDindex-IDindex2),FitInd,...
                           GA.Progress)),1,'first');      
            TopTopIDs(IDindex2) = ThisID;
            % Remove it from TopIDs
            TopIDs(TopIDs==ThisID) = [];

            IDindex2 = IDindex2+1;
            FitIndex = FitIndex+1;
            if FitIndex > GA.NFit
                FitIndex = Fit0;
            end
        end
    end
        
end

