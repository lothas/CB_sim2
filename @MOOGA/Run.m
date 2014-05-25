function [ GA ] = Run( GA )
%RUN Runs the optimization algorithm
%   The Run method creates a random population (or loads an existing
%   one), and then goes through the optimization process until the
%   last generation.

if GA.Sim.Graphics == 1
    GA.Sim.Fig = figure();
end
        
Sim = GA.Sim;
Gen = GA.Gen;
NFit = GA.NFit;
FitFcn = GA.FitFcn;
for g = GA.Progress:GA.Generations    
    gSeqs = GA.Seqs(:,:,g);
    gFit = GA.Fit(:,:,g);
    for i = 1:GA.Population
        if any(gFit(i,:)~=0)
            continue;
        end
        
        % Set-up the simulation
        wSim = copy(Sim);
        wSim = Gen.Decode(wSim,gSeqs(i,:)); %#ok<PFBNS>
        wSim = wSim.Init();
        
        % Run the simulation
        cla % clear previous render
        wSim = wSim.Run();
        
        % Calculate the genome's fitness
        thisFit = zeros(1,NFit);
        for f = 1:NFit
            thisFit(f) = FitFcn{f}(wSim); %#ok<PFBNS>
            
            if strcmp(func2str(FitFcn{f}),...
                '@(varargin)GA.VelFit(varargin{:})')
                % Switch direction if the model walks backwards
                if thisFit(f)<0
                    gSeqs(i,:) = Gen.SwitchDir(gSeqs(i,:));
                end                    
            end
        end        
        gFit(i,:) = thisFit;
    end
    GA.Fit(:,:,g) = gFit;
    GA.Seqs(:,:,g) = gSeqs;
    
    % Finished running tests, create new generation
    TopIDs = GA.GetTopPop(GA.Fittest(1)); % fitness = genes
    
    % Transfer top IDs to new population
    GA.Seqs(1:GA.Fittest(1),:,GA.Progress+1) = ...
        GA.Seqs(TopIDs,:,GA.Progress);
    
    % Add a mutated copy of top IDs
    GA.Seqs(GA.Fittest(1)+1:GA.Fittest(1)+GA.Fittest(2),:,GA.Progress+1) = ...
        GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,GA.Progress+1));
    
    % Add children (by pairs)
    for c = 1:ceil(GA.Fittest(3)/2)
        % Select 2 parents randomly
        IDs = randsample(1:GA.Fittest(1),2);
        Child = GA.Fittest(1)+GA.Fittest(2) + 2*c-1;
        if Child == GA.Population
            Children = ...
                GA.Gen.Crossover(GA.Seqs(IDs(1),:,GA.Progress+1),...
                                 GA.Seqs(IDs(2),:,GA.Progress+1));
            GA.Seqs(Child,:,GA.Progress+1) = Children(1,:);
        else
            GA.Seqs(Child:Child+1,:,GA.Progress+1) = ...
                GA.Gen.Crossover(GA.Seqs(IDs(1),:,GA.Progress+1),...
                                 GA.Seqs(IDs(2),:,GA.Progress+1));
        end
    end
    
    % Finished processing generation
    GA.Progress = g+1;
    
    % Display top results
    TopFit = GA.Fit(TopIDs,:,g);
    disp(max(TopFit));
    
    if ~isempty(GA.FileOut)
        % Save progress to file
        save(GA.FileOut,'GA');
    end
    
    if isa(GA.GenerationFcn,'function_handle')
        % Call external function
        GA.GenerationFcn();
    end
end

end

