function [ TopIDs, Weights ] = GetTopPop( GA, N )
%GETTOPPOP Receives a population's fitness and returns the N top IDs
%   GetTopPop divides the given genomes into Pareto fronts based
%   on their fitness. It then selects the top N genomes from the
%   layers one by one. If the layer is larger than required,
%   genomes that excell at each fitness individually are picked.

    if GA.NFit == 1
        FitID = [GA.Fit(:,:,GA.Progress), (1:GA.Population)'];
        Sorted = sortrows(FitID,-1);
        TopIDs = Sorted(1:N,2);
        Weights = ones(N,1);
    else
        % Divide population into pareto fronts
        Fronts = GA.Pareto(GA.Fit(:,:,GA.Progress));
        f = 1;
        IDindex = 1;
        TopIDs = zeros(GA.Population,1);
        Weights = zeros(GA.Population,1);
        while IDindex<=N
            % Take IDs from pareto fronts until TopPop is reached
            res = length(Fronts{f});
            TopIDs(IDindex:IDindex+res-1) = Fronts{f}(randperm(res));
            Weights(IDindex:IDindex+res-1) = repmat(f,res,1);
            % Fronts are randomly permutated to remove biases
            IDindex = IDindex+res;
            f = f+1;
        end

        if f==2 && IDindex>N
            % The first pareto front had more than N individuals
            % Select the best genomes from each fitness
    %         if min(Fit(1,TopIDs(1:IDindex-1)))>0.95
    %             % Don't select genomes by height fitness (Fit(1))
    %             Fit0 = 2;
    %         else
                Fit0 = 1;
    %         end

            IDindex2 = 1;
            FitIndex = 1;
            TopTopIDs = zeros(N,1);
            while IDindex2<=N
                % Select the highest genomes for each fitness
                ThisID = find(GA.Fit(:,FitIndex,GA.Progress) == ...
                    max(GA.Fit(TopIDs(1:IDindex-IDindex2),FitIndex,...
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

            TopIDs = TopTopIDs;
        end

        TopIDs = TopIDs(1:N);
        % Invert the weights (first fronts get higher weight)
        Weights = f-Weights;
        Weights = Weights(1:N);
    end
end

