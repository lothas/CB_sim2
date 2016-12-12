function [topFits, topSeqs, filenames] = Our_GA_Results()
%OUR_GA_RESULTS Summary of this function goes here
%   Detailed explanation goes here
% Find all results files
res_files = dir('VGAO*.mat');
nFiles = length(res_files);

% All GAs are assumes to have the same number of genes and fitness values
nFits = 11; nGenes = 22;
VelRangeCol = 3; % Column ID containing the VelRange fitness
nTopPerFile = 5;

topData = zeros(nFiles*nTopPerFile, nFits+nGenes+1);

for i = 1:nFiles
    % Load the results file
    load(res_files(i).name);
    
    % Show properties
    disp(['***-------------------------- ', ...
        res_files(i).name,' --------------------------***']);
    disp(['Population: ',int2str(GA.Population),...
        '   |   Generations: ',int2str(GA.Generations),...
        '   |   Runtime: ',GA.PrintRuntime()]);
    % Show fitness
    fitStr = 'Fitness functions considered: ';
    for f = 1:GA.NFit
        fitStr = [fitStr, ...
            GA.GetFitFcnName(GA.FitFcn{f,2}), ', ']; %#ok<AGROW>
    end
    disp(fitStr(1:end-2));
    GA.Find();
    
    % Get best VelRange fit genomes
    genData = [GA.Fit(:,:,GA.Progress), ...
               GA.Seqs(:,:,GA.Progress)];
	genData = sortrows(genData,-VelRangeCol);
    % Print best
    disp('Best VelRange fit genomes:');
    disp(genData(1:nTopPerFile,GA.FitIDs));
    topData((i-1)*nTopPerFile+1:i*nTopPerFile,:) = ...
        [genData(1:nTopPerFile,:), i*ones(nTopPerFile,1)];
    
    disp(['-----------------------------', ...
        '----------------------','-----------------------------']);
    disp(' ');
end

% Sort out data by best VelRange fit
topData = sortrows(topData,-VelRangeCol);
topFits = topData(:,1:nFits);
topSeqs = topData(:,nFits+1:end-1);
filesIDs = topData(:,end);
% Transform file IDs into filenames
filenames = cell(size(filesIDs));
for i = 1:length(filenames)
    filenames{i} = res_files(filesIDs(i)).name;
end

end