function [ GA ] = Run( GA )
%RUN Runs the optimization algorithm
%   The Run method creates a random population (or loads an existing
%   one), and then goes through the optimization process until the
%   last generation.

if GA.Sim.Graphics == 1
    GA.Sim.Fig = figure();
end

if isempty(gcp('nocreate'))
    parpool(7) % Work in parallel to finish faster
end

% Decode base genome if provided
if ~isempty(GA.BaseSeq)
    GA.Sim = GA.BaseGen.Decode(GA.Sim, GA.BaseSeq);
end

Sim = GA.Sim;
Gen = GA.Gen;
NFit = GA.NFit;
FitFcn = GA.FitFcn;
lastTime = [0 0];
for g = GA.Progress+1:GA.Generations
    startTime = now();
    if all(lastTime == 0)
        disp(['Running generation ',int2str(g),' of ',int2str(GA.Generations),...
        '   -   Start: ',datestr(startTime,'HH:MM')]);
    else
        if lastTime(2) == 0
            ETA = addtodate(startTime,ceil(lastTime(1)),'second');
        else
            ETA = addtodate(startTime,ceil(2*lastTime(1)-lastTime(2)),'second');
        end
        disp(['Running generation ',int2str(g),' of ',int2str(GA.Generations),...
        '   -   Start: ',datestr(startTime,'HH:MM'),' - ETA: ',datestr(ETA,'HH:MM')]);
    end    
    inner_tic = tic;
    
    gSeqs = GA.Seqs(:,:,g);
    gFit = GA.Fit(:,:,g);
    parfor i = 1:GA.Population
%     for i = 1:GA.Population
        if any(gFit(i,:)~=0)
            continue;
        end
        
        % Set-up the simulation
        wSim = deepcopy(Sim);
        wSim = Gen.Decode(wSim,gSeqs(i,:)); %#ok<PFBNS>
        wSim = wSim.Init();
        
        % Run the simulation
        wSim = wSim.Run();
        
        % Calculate the genome's fitness
        thisFit = zeros(1,NFit);
        thisOuts = cell(1,NFit);
        for f = 1:NFit
            % Preprocessing for ZMPFit
            if ~isempty(strfind(func2str(FitFcn{f}),'ZMPFit')) %#ok<PFBNS>
                % Prepare all the required vectors
                % (torques, state, etc) and put them in wSim.Out
                
                % ZMP Fit should be the last one as it uses the
                % output from all previous fitness functions
                Outs = find(~cellfun(@isempty,thisOuts));
                for o = 1:length(Outs)
                    wSim.Out = wSim.JoinOuts(thisOuts{Outs(o)});
                end
                wSim.Con.FBType = 2;
            end
            
            % Call the fitness function
            [thisFit(f),thisOuts{f}] = FitFcn{f}(wSim);
            
            % Postprocessing for VelFit
            if ~isempty(strfind(func2str(FitFcn{f}),'VelFit'))
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
            
            % Postprocesing for ZMPFit
            if ~isempty(strfind(func2str(FitFcn{f}),'ZMPFit'))
                % Update the fitness based on the uphill (4) and
                % downhill (5) fitness
                thisFit(f) = 7.5*thisFit(f)*...
                    (1-cosd(max(thisFit(4),thisFit(5))));
            end
        end        
        gFit(i,:) = thisFit;
    end
    GA.Fit(:,:,g) = gFit;
    GA.Seqs(:,:,g) = gSeqs;
    
    % Display top results
    disp(['Generation ',num2str(g),' results: ',num2str(max(GA.Fit(:,:,g)))]);
    
    % Display time elapsed
    t_diff = toc(inner_tic);
    minutes = floor(t_diff/60);
    seconds = mod(t_diff,60);
    if minutes>0
        disp(['Time to run this generation: ',num2str(minutes),''' ',num2str(seconds,'%.2f'),'"',10])
    else
        disp(['Time to run this generation: ',num2str(seconds,'%.2f'),'"',10])
    end
    lastTime(2) = lastTime(1);
    lastTime(1) = t_diff;
    
    % Finished processing generation
    GA.Progress = g;
    
    if g<GA.Generations
        % Finished running tests, create new generation
        [TopIDs,Weights] = GA.GetTopPop(GA.Fittest(1)); % fitness = genes

        % Transfer top IDs to new population
        GA.Seqs(1:GA.Fittest(1),:,g+1) = ...
            GA.Seqs(TopIDs,:,g);
        GA.Fit(1:GA.Fittest(1),:,g+1) = ...
            GA.Fit(TopIDs,:,g);

        % Add a mutated copy of top IDs
        GA.Seqs(GA.Fittest(1)+1:GA.Fittest(1)+GA.Fittest(2),:,g+1) = ...
            GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,g+1));

        % Add children (by pairs)
        Child = GA.Fittest(1)+GA.Fittest(2)+1;
        while Child<=GA.Population-GA.Fittest(3)
            % Select 2 parents randomly
            if GA.WeightedPairing
                IDs = randsample(1:GA.Fittest(1),2,'true',Weights);
                if IDs(1) == IDs(2)
                    continue
                end
            else
                IDs = randsample(1:GA.Fittest(1),2);
            end
            
            [Children, Res] = ...
                GA.Gen.Crossover(GA.Seqs(IDs(1),:,g+1),...
                                 GA.Seqs(IDs(2),:,g+1));
            for c = 1:2
                if Res{c,1} ~= 0 && Child <= GA.Population
                    GA.Seqs(Child,:,g+1) = Children(c,:);
                    Child = Child+1;
                end
            end
        end
        
        % Add random children
        GA.Seqs(GA.Population-GA.Fittest(3)+1:GA.Population,:,g+1) = ...
           GA.Gen.RandSeq(GA.Fittest(3)); 
    end
        
    if ~isempty(GA.FileOut)
        % Save progress to file
        save(GA.FileOut,'GA');
    end
    
    if isa(GA.GenerationFcn,'function_handle')
        % Call external function
        GA = GA.GenerationFcn(GA);
    end
end

end

