function [ GA ] = Run( GA )
%RUN Runs the optimization algorithm
%   The Run method creates a random population (or loads an existing
%   one), and then goes through the optimization process until the
%   last generation.

if GA.Sim.Graphics == 1
    GA.Sim.Fig = figure();
end

if isempty(gcp('nocreate'))
    parpool(6) % Work in parallel to finish faster
end

% Decode base genome if provided
if ~isempty(GA.BaseSeq)
    GA.Sim = GA.BaseGen.Decode(GA.Sim, GA.BaseSeq);
end

% Sim = GA.Sim;
% Gen = GA.Gen;
% NFit = size(GA.FitFcn,1);
% FitInd = GA.FitFcn(:,1); % Index for fit functions
% FitFcn = GA.FitFcn(:,2); % Handles for fit functions
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
    
	gMRT = GA.MLseqRunTime(:,:,g);	% Matsuoka sim runTime
	gSRT = GA.simRunTime(:,:,g);	% Matsuoka + CB runTime
    gSeqs = GA.Seqs(:,:,g);
    gFits = GA.Fit(:,:,g);
    sim_endCond = GA.sim_endCond(:,:,g);
    Tend_ratio = GA.Tend_ratio(:,:,g);
    
    ParRunSeq = @GA.RunSeq;
    
    parfor i = 1:GA.Population
%     for i = 1:GA.Population
        if any(gFits(i,:)~=0)
            continue;
        end
                
        [gFits(i,:), gSeqs(i,:), gMRT(i,1),...
            gSRT(i,1),sim_endCond(i,1),Tend_ratio(i,1)] =...
            feval(ParRunSeq, gSeqs(i,:));
        

    end
	
	GA.MLseqRunTime(:,:,g) = gMRT;
	GA.simRunTime(:,:,g) = gSRT;
    GA.Fit(:,:,g) = gFits;
    GA.Seqs(:,:,g) = gSeqs;
    GA.sim_endCond(:,:,g) = sim_endCond;
    GA.Tend_ratio(:,:,g) = Tend_ratio;
    
    % Check outliers
    ids = [];
    for n = 1:length(GA.FitIDs)
        ids = [ids; find(GA.Fit(:,n,g) > ...
                    mean(GA.Fit(:,n,g))+3.5*std(GA.Fit(:,n,g)))]; %#ok<AGROW>
    end
    ids = unique(ids);
    newFits = GA.Fit(ids,:,g);
    parfor i = 1:length(ids)  
        [newFits(i,:), ~] = feval(ParRunSeq, gSeqs(i,:));
    end
    GA.Fit(ids,:,g) = min(GA.Fit(ids,:,g), newFits);
    
    % Display top results
%     PFits = repmat(GA.FitMinMax, size(GA.Fit(:,:,g), 1), 1).*GA.Fit(:,:,g);
%     MaxFits = GA.FitMinMax.*max(PFits);
%     disp(['Generation ',num2str(g),' results: ',num2str(MaxFits, '  %.4f')]);
    disp(['Generation ',num2str(g),' results:']);
    GA.Find();
    
    % Display time elapsed
    t_diff = toc(inner_tic);
    GA.totGenTime(1,g) = t_diff; % save the generation time.
    
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
    GA.CompTime = GA.CompTime + t_diff;
    
    if g<GA.Generations
        % Finished running tests, create new generation
        [TopIDs,Weights] = GA.GetTopPop(GA.Fittest(1)); % fitness = genes

        % Transfer top IDs to new population
        GA.Seqs(1:GA.Fittest(1),:,g+1) = ...
            GA.Seqs(TopIDs,:,g);
        GA.Fit(1:GA.Fittest(1),:,g+1) = ...
            GA.Fit(TopIDs,:,g);
        GA.Parents(1:GA.Fittest(1),:,g+1) = ...
            GA.Parents(TopIDs,:,g);

        % Add a mutated copy of top IDs
        GA.Seqs(GA.Fittest(1)+1:GA.Fittest(1)+GA.Fittest(2),:,g+1) = ...
            GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,g+1));
        GA.Parents(GA.Fittest(1)+1:GA.Fittest(1)+GA.Fittest(2),:,g+1) = ...
            GA.Parents(GA.Fittest(1)+1:GA.Fittest(1)+GA.Fittest(2),:,g);

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
                    GA.Parents(Child,:,g+1) = IDs;
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

