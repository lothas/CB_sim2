function MOGAcrowd(mode, type)
% GATEST Runs a MOGA test-case.
%   A MOGA is run using the following genome:
%   [x, y] - point position in (x,y) coordinates
%   and the following fitness functions:
%   Dist_x: distance of x from 0
%   Dist_y: distance of y from 0
%   the genome is bounded by x^2+y^2=1
%   The algorithm has different selection modes:
%   1 - weighted selection (equal weights)
%   2 - max genomes are selected from each fitness
%   3 - pareto front selection
%   4 - pareto front + crowding selection
%   Type: 1 - square, 2 - circle
if nargin<2
    type = 2;
    if nargin<1
        mode = 5;
    end
end

% Algorithm parameters
Pop = 200;
NGen = 100;
NTop = 40;

AllPop = zeros(Pop,2,NGen);
AllTop = zeros(NTop,2,NGen);
AllFit = zeros(Pop,2,NGen);
AllDistMu = zeros(NGen,1);
AllDistMin = zeros(NGen,1);
stdist = zeros(NGen,1);

% Create random population
GenPop = RandomG(Pop);

percExtra = 1.1;
curGen = 0;
while curGen<NGen
    % Obtain population fitness
    Fit = zeros(Pop,2);
    Fit(:,1) = Fit1(GenPop);
    Fit(:,2) = Fit2(GenPop);
    weights = ones(Pop,1)/Pop;
    
    % Store results
    AllPop(:,:,curGen+1) = GenPop;
    AllFit(:,:,curGen+1) = Fit;
    
    if mode == 1
        % Use weighted selection method
        wFit = Fit(:,1)+Fit(:,2);
        wFitID = [wFit,(1:Pop)'];
        wFitID = sortrows(wFitID,-1);
        BestIDs = wFitID(1:NTop,end);
    else
        if mode == 2
            % Use max selection method
            FitID = [Fit,(1:Pop)'];
            FitID1 = sortrows(FitID,-1);
            FitID2 = sortrows(FitID,-2);
            BestIDs = [FitID1(:,end),FitID2(:,end)];
            BestIDs = reshape(BestIDs',[],1);
            BestIDs = BestIDs(1:NTop);
        else
            % Sort into Pareto fronts
            Fronts = Pareto(Fit);
            
            if mode == 3
                % Select from Pareto fronts until NTop
                f = 1;
                BestIDs = [];
                while length(BestIDs)<NTop
                    BestIDs = [BestIDs;Fronts{f}]; %#ok<AGROW>
                    f = f + 1;
                end
                BestIDs = BestIDs(1:NTop);
                
                % Update distance between all candidates
                CFits = Fit(BestIDs,:);
                NC = size(CFits,1);
                dists = zeros(NC,1);
                for c = 1:NC
                    distM = repmat(CFits(c,:),NC-1,1) - ...
                        CFits([1:c-1,c+1:end],:);
                    dist = diag(distM*distM');
                    dists(c) = min(dist);
                end
                weights = dists;
            end
            
            if mode == 4
                % Mode 4 selects the best genomes based on
                % pareto level and relative distance
                f = 1;
                BestIDs = [];
                Candidates = [];
                while length(BestIDs)<NTop
                    % Add next front to Candidates
                    Candidates = [Candidates;
                                  Fronts{f},zeros(length(Fronts{f}),1)];
                    f = f + 1;
                    CanPop = Fit(Candidates(:,1),:);
                    
                    % Check distance between all candidates
                    NC = size(Candidates,1);
                    if NC > 2
                        dists = zeros(NC,1);
                        for c = 1:NC
                            distM = repmat(CanPop(c,:),NC-1,1) - ...
                                CanPop([1:c-1,c+1:end],:);
                            dist = diag(distM*distM');
                            dists(c) = min(dist);
                        end

                        % Move candidates with high variety to BestIDs
                        cIDs = dists>mean(dists)-curGen/NGen*std(dists);
%                         cIDs = dists>mean(dists)-std(dists);
%                         cIDs = dists>mean(dists);
%                         BestC = Candidates(cIDs);
%                         Candidates(cIDs) = [];
%                         BestIDs = [BestIDs;BestC];
                        Candidates(cIDs,2) = 1;
                        if sum(Candidates(:,2))>=NTop
                            BestIDs = Candidates(Candidates(:,2)==1,1);
                        end
                    end
                end
                BestIDs = BestIDs(1:NTop);
                
                % Update distance between all candidates
                CFits = Fit(BestIDs,:);
                NC = size(CFits,1);
                dists = zeros(NC,1);
                for c = 1:NC
                    distM = repmat(CFits(c,:),NC-1,1) - ...
                        CFits([1:c-1,c+1:end],:);
                    dist = diag(distM*distM');
                    dists(c) = min(dist);
                end
                weights = dists;
            end
            
            if mode == 5
                % Select from Pareto fronts until NTop+10%
                f = 1;
                BestIDs = [];
                while length(BestIDs)<NTop
                    BestIDs = [BestIDs;Fronts{f}]; %#ok<AGROW>
                    f = f + 1;
                end
                
                % Check distance between all candidates
                CFits = Fit(BestIDs,:);
                NC = size(CFits,1);
                dists = zeros(NC,1);
                for c = 1:NC
                    distM = repmat(CFits(c,:),NC-1,1) - ...
                        CFits([1:c-1,c+1:end],:);
                    dist = diag(distM*distM');
                    dists(c) = min(dist);
                end
                stdist(curGen+1) = std(dists)/mean(dists);
                dists = [dists,BestIDs];
                
                while length(BestIDs)>NTop
                    % Find pair of closest points
                    mIDs = find(dists(:,1) == min(dists(:,1)));
                    dIDs = mIDs;
                    if length(dIDs)>=2
                        % Keep only one point at random
                        safeID = randsample(length(dIDs),1);
                        dIDs(safeID) = [];
                    end
                    BestIDs(ismember(BestIDs,dists(dIDs,2))) = [];
                    dists(mIDs,:) = [];
                end
%                 BestIDs = BestIDs(1:NTop);
                
                % Update distance between all candidates
                CFits = Fit(BestIDs,:);
                NC = size(CFits,1);
                dists = zeros(NC,1);
                for c = 1:NC
                    distM = repmat(CFits(c,:),NC-1,1) - ...
                        CFits([1:c-1,c+1:end],:);
                    dist = diag(distM*distM');
                    dists(c) = min(dist);
                end
                weights = dists;
            end
        end
    end
    
    % Store results
    TopFit = Fit(BestIDs,:);
    AllTop(:,:,curGen+1) = TopFit;
    dists = zeros(NTop,1);
    for c = 1:NTop
        distM = repmat(TopFit(c,:),NTop-1,1) - ...
            TopFit([1:c-1,c+1:end],:);
        dist = diag(distM*distM');
        try
            dists(c) = min(dist);
        catch
            disp(1);
        end
    end
    AllDistMu(curGen+1) = mean(dists);
    AllDistMin(curGen+1) = min(dists);
    
    % Build new population
    NewPop = zeros(size(GenPop));
    NewPop(1:NTop,:) = GenPop(BestIDs,:); % Elitism
    g = NTop+1;
    while g<Pop
        Parents = randsample(BestIDs,2,true,weights);
        if Parents(1)~=Parents(2)
            go = Offspring(GenPop(Parents(1),:),...
                           GenPop(Parents(2),:));
            if ~isempty(go)
                NewPop(g,:) = go;
                g = g + 1;
            end
        end
    end
    
    GenPop = NewPop;
    curGen = curGen + 1;
%     disp (curGen);
end

% Display results
figure
x = cos(0:0.01:pi/2);
y = sin(0:0.01:pi/2);
plot(x,y,'--r','LineWidth',2);
hold on
for g = 1:NGen
    Color = ((NGen-g)*[1, 0, 0]+g*[0.1, 0.8, 0.2])/NGen;
    scatter(AllTop(:,1,g),AllTop(:,2,g),50,Color,'fill');
%     scatter(AllFit(:,1,g),AllFit(:,2,g),50,Color,'fill');
end
axis equal
axis([0 1 0 1])

figure
plot(stdist);
% figure
% plot(AllDistMu);
% hold on
% plot(AllDistMin);
% disp([AllDistMu(end), AllDistMin(end)]);

    function g = RandomG(n)
        radii = rand(n,1);
        angle = pi/2*rand(n,1);
        g = [radii.*cos(angle),radii.*sin(angle)];
    end

    function f = Fit1(g)
        f = g(:,1);
    end

    function f = Fit2(g)
        f = g(:,2);
    end

    function g = Offspring(p1,p2)
        % Crossover
        if rand()>0.5
            g = [p1(1),p2(2)];
        else
            g = [p2(1),p1(2)];
        end
        
        % Mutation
        g = g + 0.2*(0.5-rand(1,2));
        
        % Bounding
%         g = min(max(g,[0, 0]),[1, 1]);
        if type == 1
            if any(g<0) || any(g>1)
                g = [];
            end
        else
            if g*g'>1
                g = [];
            end
        end
    end

    function Fronts = Pareto(Data)
        Npoints = size(Data,1);
        
        % Add IDs to data points
        X = [Data,(1:Npoints)'];
        
        Fronts = {};
        fn = 0;
        XTemp = [X,zeros(size(X,1),1)];
        Under = [];
        while ~isempty(XTemp)
            p = find(XTemp(:,end)==0,1,'first');
            if isempty(p)
                % Finished checking the front
                fn = fn + 1;
                Front = XTemp(:,end-1);
                Fronts{fn,1} = Front;
                
                % Prepare data to evaluate next front
                XTemp = [Under,zeros(size(Under,1),1)];
                Under = [];
            else
                XTemp(p,end) = 1;  % set p as "checked"
                % Find which items in XTemp are dominated by p
                pX = repmat(XTemp(p,1:end-2),size(XTemp,1),1);
                IDs = any(pX>XTemp(:,1:end-2),2) & ...
                    all(pX>=XTemp(:,1:end-2),2);
                % Move those to Under
                Under = [Under;XTemp(IDs,1:end-1)];
                XTemp(IDs,:) = [];
            end
        end     
    end
end

