function [MutSeq,Res] = Mutate(Ge,Seq,Try)
%MUTATE Creates mutations of the sequences provided
%   Mutate takes a number of sequences and adds mutations
%   as specified by the Genome object parameters:
%   MutProb: single gene mutation probability
%   MutDelta: Max strength of mutation as % of allowed range
%   MutType: Distribute mutation strength as uniform or normal distribution
    if nargin<3
        Try = 1;
    end

    [NS,NG] = size(Seq);
    if NG~=Ge.Length
        error(['Decode failed. Genes provided: ',num2str(length(Seq)),...
            '. Genes required: ',num2str(Ge.Length)]);
    end

    % Find which genes will mutate
    doMutate = rand(NS,NG)<Ge.MutProb;
    % Set the strength of each mutation
    range = Ge.MutDelta*repmat(Ge.Range(2,:)-Ge.Range(1,:),NS,1);
    switch Ge.MutType
        case {'uni','uniform'}
            Mutations = (1-2*rand(NS,NG)).*range;
        case {'norm','normal'}
            Mutations = randn(NS,NG).*range/6;
        otherwise
            error(['Mutation type ',Ge.MutType,' unknown']);
    end
    MutSeq = Seq + doMutate.*Mutations;

    % Check the mutated sequences
    Res = cell(NS,2);
    for s = 1:NS
        [thisRes,MutSeq(s,:)] = ...
            Ge.CheckGenome(MutSeq(s,:));
        if all(MutSeq(s,:) == Seq(s,:)) && thisRes{1} ~= 0
            % The sequence didn't mutate at all
            Res(s,:) = {-1,'Sequence left unchanged'};
        else
            Res(s,:) = thisRes;
        end
        
        if Res{s,1} == 0
            switch Try
                case {1,2,3,4,5,6,7}
                    % Mutation failed, try again
                    [MutSeq(s,:),Res(s,:)] = Ge.Mutate(Seq(s,:),Try+1);                    
                case 8
                    % display the problematic sequence
                    disp(Ge.seq2str(Seq(s,:)))
                case 9
                    % Try bounding the original sequence
                    disp('Mutate bounded')
                    BoundSeq = max(min(Seq(s,:),Ge.Range(2,:)),Ge.Range(1,:));
                    [MutSeq(s,:),Res(s,:)] = Ge.Mutate(BoundSeq,Try+1);
                case 10
                    % Stop trying and return original sequence
                    MutSeq(s,:) = Seq(s,:);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [thisRes,~] = ...
                        Ge.CheckGenome(MutSeq(s,:));
                    if thisRes{1}~=1
                        disp(['Mutate: ',thisRes{2}])
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end