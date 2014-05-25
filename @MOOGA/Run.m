function [ GA ] = Run( GA )
%RUN Runs the optimization algorithm
%   The Run method creates a random population (or loads an existing
%   one), and then goes through the optimization process until the
%   last generation.

if GA.Sim.Graphics == 1
    GA.Sim.Fig = figure();
end

if matlabpool('size')==0
    matlabpool open 7 % Work in parallel to finish faster
end
    
Sim = GA.Sim;
Gen = GA.Gen;
NFit = GA.NFit;
FitFcn = GA.FitFcn;
for g = GA.Progress:GA.Generations    
    gSeqs = GA.Seqs(:,:,g);
    gFit = GA.Fit(:,:,g);
    parfor i = 1:GA.Population
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
                    revSeq = Gen.SwitchDir(gSeqs(i,:));
                    [Res,revSeq] = Gen.CheckGenome(revSeq);
                    if Res{1}
                        gSeqs(i,:) = revSeq;
                        thisFit(f) = -thisFit(f);
                    end
                end                    
            end
        end        
        gFit(i,:) = thisFit;
    end
    GA.Fit(:,:,g) = gFit;
    GA.Seqs(:,:,g) = gSeqs;
    
    if g<GA.Generations
        % Finished running tests, create new generation
        TopIDs = GA.GetTopPop(GA.Fittest(1)); % fitness = genes

        % Transfer top IDs to new population
        GA.Seqs(1:GA.Fittest(1),:,GA.Progress+1) = ...
            GA.Seqs(TopIDs,:,GA.Progress);

        % Add a mutated copy of top IDs
        GA.Seqs(GA.Fittest(1)+1:GA.Fittest(1)+GA.Fittest(2),:,GA.Progress+1) = ...
            GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,GA.Progress+1));

        % Add children (by pairs)
        Child = GA.Fittest(1)+GA.Fittest(2)+1;
        while Child<=GA.Population
            % Select 2 parents randomly
            IDs = randsample(1:GA.Fittest(1),2);
            [Children, Res] = ...
                GA.Gen.Crossover(GA.Seqs(IDs(1),:,GA.Progress+1),...
                                 GA.Seqs(IDs(2),:,GA.Progress+1));
            for c = 1:2
                if Res{c,1} ~= 0 && Child <= GA.Population
                    GA.Seqs(Child,:,GA.Progress+1) = Children(c,:);
                    Child = Child+1;
                end
            end
        end

        % Finished processing generation
        GA.Progress = g+1;
    end
    
    % Display top results
    TopFit = GA.Fit(TopIDs,:,g);
    disp(['Generation ',num2str(g),' results: ',num2str(max(TopFit))]);
    
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

