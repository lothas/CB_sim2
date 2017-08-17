function [ GA ] = InitGen( GA )
%INIT Initializes the optimization variables
%   If an input MOOGA isn't provided, a random first generation is created

% Initialize sequences and fitness
GA.Seqs = zeros(GA.Population, GA.Gen.Length, GA.Generations);
GA.Parents = zeros(GA.Population, 2, GA.Generations); % 2 parents
GA.Fit = zeros(GA.Population, max(cell2mat(GA.FitFcn(:,1)')), GA.Generations);
GA.MLseqRunTime = NaN(GA.Population, 1, GA.Generations);
GA.simRunTime = NaN(GA.Population, 1, GA.Generations);
GA.sim_endCond = zeros(GA.Population, 1, GA.Generations);
GA.Tend_ratio = zeros(GA.Population, 1, GA.Generations);

GA.totGenTime = zeros(1,GA.Generations);
% TODO: do for 'MLseqRunTime' and 'simRunTime' what is done for 'Fit' and
%       'Seqs' (' if exist('In','var') == 1 ')

if exist(GA.FileIn,'file') == 2
    % Load input file
    In = load(GA.FileIn);
else
    if exist(GA.FileOut,'file') == 2
        % Load input file
        In = load(GA.FileOut);
    end
end

if exist('In','var') == 1
    % Are the same fitness functions used?
    InNames = cell(In.GA.NFit,1);
    ThisFitID = [];
    InFitID = [];
    for i = 1:In.GA.NFit
        InNames{i} = MOOGA.GetFitFcnName(In.GA.FitFcn{i,2});
    end
    for i = 1:GA.NFit
        ThisName = MOOGA.GetFitFcnName(GA.FitFcn{i,2});
        ID = find(strcmp(ThisName,InNames),1,'first');
        if ~isempty(ID)
            ThisFitID = [ThisFitID,GA.FitFcn{i,1}]; %#ok<AGROW>
            InFitID = [InFitID,In.GA.FitFcn{ID,1}]; %#ok<AGROW>
        end
        
        if length(ThisFitID) ~= length(InFitID)
            error('Discrepancy in FitFcn indexes')
        end
    end
    
    if In.GA.Population == GA.Population
        if GA.ReDo
            % Copy the last generation's seq. into new GA
            % ids = [1,3:19];
            % GA.Seqs(:,:,1) = In.GA.Seqs(:,ids,In.GA.Progress);
            GA.Seqs(:,:,1) = In.GA.Seqs(:,:,In.GA.Progress);
            GA.Fit(:,ThisFitID,1) = ...
                In.GA.Fit(:,InFitID,In.GA.Progress);
            GA.Parents(:,:,1) = In.GA.Parents(:,:,In.GA.Progress);
            GA.MLseqRunTime(:,1,1) = ...
                In.GA.MLseqRunTime(:,1,In.GA.Progress);
            GA.simRunTime(:,1,1) = ...
                In.GA.simRunTime(:,1,In.GA.Progress);
            GA.sim_endCond(:,1,1) = ...
                In.GA.sim_endCond(:,1,In.GA.Progress);
            GA.Tend_ratio(:,1,1) = ...
                In.GA.Tend_ratio(:,1,In.GA.Progress);
            GA.totGenTime(1,1) = ...
                In.GA.totGenTime(1,In.GA.Progress);
        else
            % Copy all the progress into new GA
            GA.Seqs(:,:,1:In.GA.Progress) = ...
                In.GA.Seqs(:,:,1:In.GA.Progress);
            GA.Fit(:,ThisFitID,1:In.GA.Progress) = ...
                In.GA.Fit(:,InFitID,1:In.GA.Progress);
            GA.Parents(:,:,1) = ...
                In.GA.Parents(:,:,1:In.GA.Progress);
            GA.MLseqRunTime(:,1,1) = ...
                In.GA.MLseqRunTime(:,1,1:In.GA.Progress);
            GA.simRunTime(:,1,1) = ...
                In.GA.simRunTime(:,1,1:In.GA.Progress);
            GA.sim_endCond(:,1,1) = ...
                In.GA.sim_endCond(:,1,1:In.GA.Progress);
            GA.Tend_ratio(:,1,1) = ...
                In.GA.Tend_ratio(:,1,1:In.GA.Progress);
            GA.totGenTime(1,1) = ...
                In.GA.totGenTime(1,1:In.GA.Progress);
            
            GA.Progress = In.GA.Progress-1;
        end
    else
        if In.GA.Population > GA.Population
            warning('this code is not changed yet to do use bigger Pop!');
            erorr('change it!');
            
            % Select the best genomes from the last generation
            In.GA.Quant = GA.Quant;
            TopIDs = In.GA.GetTopPop(GA.Population); % fitness = genes

            % Transfer top IDs to new population
            GA.Seqs(:,:,1) = In.GA.Seqs(TopIDs,:,In.GA.Progress);  
            if ~GA.ReDo 
                GA.Fit(:,ThisFitID,1) = ...
                    In.GA.Fit(TopIDs,InFitID,In.GA.Progress);
            end
        else
            % Copy the last generation's seq. into new GA
            GA.Seqs(1:In.GA.Population,:,1) = ...
                In.GA.Seqs(:,:,In.GA.Progress);
            if ~GA.ReDo
                GA.Fit(1:In.GA.Population,ThisFitID,1) = ...
                    In.GA.Fit(:,InFitID,In.GA.Progress);
            end
            % Generate new random sequences
            GA.Seqs(In.GA.Population+1:end,:,1) = ...
                GA.Gen.RandSeq(GA.Population-In.GA.Population);
        end
    end
else
    % Generate new random sequences
    GA.Seqs(:,:,1) = GA.Gen.RandSeq(GA.Population);
%     GA.Seqs(1,:,1) = [1.378073592365504   0.601677928713766  -0.153495857899771,...
%                       0.005000000000000  -1.997955184234347   0.246232104453995,...
%                       0.520031614524507  11.201294534323361   0.012676038166924   0.035076340435940];
%     GA.Seqs(1,:,1) = [1.307792274230613   0.618822719574968  -0.120891611956398,...
%         0.005000000000000  -1.925739666296446   0.256751823628962,...
%         0.509645655286311  10.503430925243782   0.045494476790015   0.066446845197164];
end
                        
end

