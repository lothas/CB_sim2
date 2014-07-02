function [ TopSeqs, TopFits ] = RescueFit( GA,N )
%RESCUEFIT 
%   Parses the fitness of all generations looking for the N
%   fittest individuals that got lost for some reason

% Reshape the matrices into one big generation
AllFit = reshape(permute(GA.Fit,[2 1 3]),GA.NFit,[],1)';
AllSeq = reshape(permute(GA.Seqs,[2 1 3]),GA.Gen.Length,[],1)';

CodedFit = cell(1,GA.NFit);
NIDs = GA.Population*GA.Generations;
for f = 1:GA.NFit
    % Prepare IDs for each genome in each fitness group
    CodedFit{f} = [(1:NIDs)', AllFit(:,f)];
    % Sort each group
    CodedFit{f} = sortrows(CodedFit{f},-2);
    % Lose the fitness
    CodedFit{f} = CodedFit{f}(:,1);
end
SortedIDs = cell2mat(CodedFit);

% Get all IDs in row and remove repeated ones
[A,~,C] = unique(reshape(SortedIDs',[],1));
TopIDs = A(C(1:N));

TopFits = AllFit(TopIDs,:);
TopSeqs = AllSeq(TopIDs,:);

end

