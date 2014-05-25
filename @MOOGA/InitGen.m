function [ GA ] = InitGen( GA )
%INIT Initializes the optimization variables
%   If an input MOOGA isn't provided, a random first generation is created

% Initialize sequences and fitness
GA.Seqs = zeros(GA.Population, GA.Gen.Length, GA.Generations);
GA.Fit = zeros(GA.Population, GA.NFit, GA.Generations);

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
    if In.GA.Population == GA.Population
        if GA.ReDo
            % Copy the last generation's seq. into new GA
            GA.Seqs(:,:,1) = In.GA.Seqs(:,:,In.GA.Progress);
        else
            % Copy all the progress into new GA
            GA.Seqs(:,:,1:In.GA.Progress) = ...
                In.GA.Seqs(:,:,1:In.GA.Progress);
            GA.Fit(:,:,1:In.GA.Progress) = ...
                In.GA.Fit(:,:,1:In.GA.Progress);
            GA.Progress = In.GA.Progress;
        end
    else
        if In.GA.Population > GA.Population
            % Select the best genomes from the last generation
            TopIDs = In.GA.GetTopPop(GA.Population); % fitness = genes

            % Transfer top IDs to new population
            GA.Seqs(:,:,1) = In.GA.Seqs(TopIDs,:,In.GA.Progress);  
            GA.Fit(:,:,1) = In.GA.Fit(TopIDs,:,In.GA.Progress);            
        else
            % Copy the last generation's seq. into new GA
            GA.Seqs(1:In.GA.Population,:,1) = ...
                In.GA.Seqs(:,:,In.GA.Progress);
            GA.Fit(1:In.GA.Population,:,1) = ...
                In.GA.Fit(:,:,In.GA.Progress);   
            % Generate new random sequences
            GA.Seqs(In.GA.Population+1:end,:,1) = ...
                GA.Gen.RandSeq(GA.Population-In.GA.Population);
        end
    end
else
    % Generate new random sequences
    GA.Seqs(:,:,1) = GA.Gen.RandSeq(GA.Population);
end

end

