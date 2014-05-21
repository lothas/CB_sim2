function [ TopIDs ] = GetTopPop( GA, Fit, N )
%GETTOPPOP Receives a population's fitness and returns the N top IDs
%   GetTopPop divides the given genomes into Pareto fronts based
%   on their fitness. It then selects the top N genomes from the
%   layers one by one. If the layer is larger than required,
%   genomes that excell at each fitness individually are picked.

    NS = size(Fit,1);
    
    % Divide population into pareto fronts
    Fronts = GA.Pareto(Fit);
    f = 1;
    IDindex = 1;
    TopIDs = zeros(NS,1);
    while IDindex<=N
        % Take IDs from pareto fronts until TopPop is reached
        res = length(Fronts{f});
        TopIDs(IDindex:IDindex+res-1) = Fronts{f}(randperm(res));
        % Fronts are randomly permutated to remove biases
        IDindex = IDindex+res;
        f = f+1;
    end
    
    if f==2 && IDindex>N
        % The first pareto front had more than N individuals
        % Select the best genomes from each fitness
        NumFits = size(Fit,1);
%         if min(Fit(1,TopIDs(1:IDindex-1)))>0.95
%             % Don't select genomes by height fitness (Fit(1))
%             Fit0 = 2;
%         else
%             Fit0 = 1;
%         end
        
        IDindex2 = 1;
        FitIndex = 1;
        TopTopIDs = zeros(N,1);
        while IDindex2<=N
            % Select the highest genomes for each fitness
            ThisID = find(Fit(FitIndex,:)==max(Fit(FitIndex,TopIDs(1:IDindex-IDindex2))),1,'first');            
            TopTopIDs(IDindex2) = ThisID;
            % Remove it from TopIDs
            TopIDs(TopIDs==ThisID) = [];
            
            IDindex2 = IDindex2+1;
            FitIndex = FitIndex+1;
            if FitIndex > NumFits
                FitIndex = Fit0;
            end
        end
        
        TopIDs = TopTopIDs;
    end
    
    TopIDs=TopIDs(1:N);
end

