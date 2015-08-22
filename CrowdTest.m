function CrowdTest()

Nsamples = 1000;
Nout = 100;

radii = rand(Nsamples,1);
angles = pi/2*rand(Nsamples,1);
X = [radii.*sin(angles), radii.*cos(angles)];
Xe = [X,(1:Nsamples)'];

Fronts = Pareto(Xe);

BestPar = [];
i = 1;
while length(BestPar)<Nout
    BestPar = [BestPar;Fronts{i}];  % best selection based on initial fit only
    i = i + 1;
end

BestCrowd = BestPar;
BestPar = BestPar(1:Nout);

% Check distance between all candidates
CFits = Xe(BestCrowd,1:end-1);
[NC,NF] = size(CFits);
CFits = [CFits,(1:NC)'];
dists = ones(NC,1);
Extr = max(CFits(:,1:end-1),[],1) - ...
    min(CFits(:,1:end-1),[],1); % Distance value given to
                                % extreme points
for f = 1:NF
    FitVal = sortrows(CFits(:,[f, end]),1);
    FitDist = diff(FitVal(:,1));
    FitDistAvg = (FitDist(1:end-1)+FitDist(2:end))/2;
    dists(FitVal(:,2)) = dists(FitVal(:,2)) .* ...
        [Extr(f); FitDistAvg; Extr(f)];
end

dists = [dists,BestCrowd];

% Remove delta candidates
while length(BestCrowd)>Nout
    % Find pair of closest points
    mIDs = find(dists(:,1) == min(dists(:,1)));
    dIDs = mIDs;
    if length(dIDs)>=2
        % Keep only one point at random
        safeID = randsample(length(dIDs),1);
        dIDs(safeID) = [];
    end
    BestCrowd(ismember(BestCrowd,dists(dIDs,2))) = [];
    dists(mIDs,:) = [];
end
                
BestPar = BestPar(1:Nout);
% % Get samples selected so far
% XF = Xe(BestPar,:);
% NCsamples = size(XF,1);
% Dist = zeros(NCsamples,1);
% for s = 1:NCsamples
%     % Calculate the distance from each sample to the others
%     FitDiff = repmat(XF(s,1:end-1),NCsamples-1,1) - ...
%         XF([1:s-1,s+1:end],1:end-1);
%     sDist = diag(FitDiff*FitDiff');
%     Dist(s) = min(sDist); % minimal distance of sample s to all others
% end
% 
% Dist = [Dist, XF(:,end)];
% while length(BestCrowd)>Nout
%     % Find the closest pair of points
%     IDs = find(Dist(:,1)==min(Dist(:,1)));
%     dIDs = IDs;
%     % Keep one point
%     safeID = randsample(length(IDs),1);
%     dIDs(safeID) = [];
%     
%     % Remove points from Dist and BestCrowd
%     BestCrowd(ismember(BestCrowd,Dist(dIDs,2))) = [];
%     Dist(IDs,:) = [];
% end

% BestCrowd = [];
% while length(BestCrowd)<Nout
%     % Get samples selected so far
%     XF = Xe(BestPar,:);
%     NCsamples = size(XF,1);
%     Dist = zeros(NCsamples,1);
%     for s = 1:NCsamples
%         % Calculate the distance from each sample to the others
%         FitDiff = repmat(XF(s,1:end-1),NCsamples-1,1) - ...
%             XF([1:s-1,s+1:end],1:end-1);
%         sDist = diag(FitDiff*FitDiff');
%         Dist(s) = min(sDist); % minimal distance of sample s to all others
%     end
%     % Select only the samples that are separated from each other by at
%     % least the mean + 1 std distance
%     PassIDs = XF(Dist>mean(Dist),end);
%     if length(PassIDs)>=Nout
%         BestCrowd = PassIDs(1:Nout);
%     else
%         i = i + 1;
%         BestPar = [BestPar;Fronts{i}];  % keep selecting samples
%     end
% end
% BestPar = BestPar(1:Nout);

% BestPar = [];
% BestCrowd = [];
% NFronts = length(Fronts);
% i = 1;
% while length(BestPar)<Nout || length(BestCrowd)<Nout
%     BestPar = [BestPar;Fronts{i}];  % best selection based on initial fit only
%     
%     % Get samples from Front i
%     XF = Xe(Fronts{i},:);
%     Dist = zeros(size(XF,1),1);
%     for s = 1:size(XF,1)
%         % Calculate the distance from each sample to the others
%         FitDiff = repmat(X(s,:),Nsamples-1,1)-X([1:s-1,s+1:end],:);
%         sDist = diag(FitDiff*FitDiff');
%         Dist(s) = min(sDist); % minimal distance of sample s to all others
%     end
%     % Select only the samples that are separated from each other by at
%     % least the mean distance
%     PassIDs = XF(Dist>mean(Dist),end);
%     BestCrowd = [BestCrowd;PassIDs];
%     i = i + 1;
% end
% BestPar = BestPar(1:Nout);
% BestCrowd = BestCrowd(1:Nout);

% BestPar = [];
% CrowdFit = zeros(size(X,1),2);
% NFronts = length(Fronts);
% for i = 1:NFronts
%     BestPar = [BestPar;Fronts{i}];  % best selection based on initial fit only
%     CrowdFit(Fronts{i},1) = NFronts - i;  % set new fit as pareto level
% end
% BestPar = BestPar(1:Nout);
% 
% 
% for i = 1:Nsamples
%     FitDiff = repmat(X(i,:),Nsamples-1,1)-X([1:i-1,i+1:end],:);
%     Dist = diag(FitDiff*FitDiff');
%     CrowdFit(i,2) = min(Dist);
% end
% 
% Fronts = Pareto([CrowdFit,(1:Nsamples)']);
% BestCrowd = [];
% NFronts = length(Fronts);
% for i = 1:NFronts
%     BestCrowd = [BestCrowd;Fronts{i}];  % best selection based on crowd
% end

figure
scatter(X(BestCrowd,1),X(BestCrowd,2),225,[0 0.8 0],'fill');
hold on
scatter(X(BestPar,1),X(BestPar,2),100,[1 0 0],'fill');
scatter(X(:,1),X(:,2),25,[0 0 0],'fill');

function Fronts = Pareto(Data,NFrontsReq)
    if nargin<2
        NFrontsReq = [];
    end
    
    % Sort by first objective
    SortedData = sortrows(Data,-1);

    if size(Data,2) == 1
        Fronts = num2cell(SortedData(:,end));
        return;
    end

    % Start separating layers
    NFronts = 1;
    Fronts = {};
    out = [];

    while 1
        % Drop samples that are "dominated"
        % i.e. that x1>x2 for each objective
        i = 1;
        Nd = size(SortedData,1);
        if Nd == 0
            break;
        end

        while i < Nd
            j = i+1;

            % Find dominated values efficiently :)
            IDs = j-1 + find(all(SortedData(j:Nd,2:end-1) <= ...
                            repmat(SortedData(i,2:end-1),Nd-j+1,1),2)==1);
            SortedData(IDs,:) = [];
            Nd = size(SortedData,1);

            i = i+1;
        end

        % Save "undominated" samples to current front
        Fronts{NFronts} = SortedData(:,end); %#ok<AGROW>
        out = [out; Fronts{NFronts}]; %#ok<AGROW>
        NFronts = NFronts+1;

        if ~isempty(NFrontsReq) && NFronts>NFrontsReq
            break;
        end

        % Restore samples to original data minus samples already on fronts
        SortedData = Data; SortedData(ismember(SortedData(:,end),out),:) = [];
        SortedData = sortrows(SortedData,-1);
    end
end

end