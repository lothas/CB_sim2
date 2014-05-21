function [MutSeq,Res] = Mutate(Ge,Seq)
%MUTATE Creates mutations of the sequences provided
%   Mutate takes a number of sequences and adds mutations
%   as specified by the Genome object parameters:
%   MutProb: single gene mutation probability
%   MutDelta: Max strength of mutation as % of allowed range
%   MutType: Distribute mutation strength as uniform or normal distribution

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
    Res = ones(NS,1);
    for s = 1:NS
        if MutSeq(s,:) == Seq(s,:)
            % The sequence didn't mutate at all
            Res(s) = -1;
        else
            [MutSeq(s,:),thisRes] = ...
                Ge.CheckGenome(MutSeq(s,:));
            Res(s) = thisRes{1};
        end

        while Res(s) == 0
            % Mutation failed, try again
            [MutSeq(s,:),Res] = Ge.Mutate(Seq(s,:));
        end
    end
end