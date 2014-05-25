function [ Seq ] = RandSeq(Ge,N,dist)
%RANDSEQ Returns N random genetic sequences
%   RandSeq creates N random genomes and returns them
%   as a (N x genome length) matrix
    switch nargin
        case 1
            N = 1;
            dist = 'uni';
        case 2
            dist = 'uni';
    end
    NGenes = size(Ge.Range,2);
    Seq = zeros(N,NGenes);
    
    for gen = 1:NGenes
        switch dist
            case {'uni','uniform'}
                Seq(:,gen) = Ge.Range(1,gen)+rand(N,1)*(Ge.Range(2,gen)-Ge.Range(1,gen));
            case {'norm','normal'}
                Seq(:,gen) = (Ge.Range(2,gen)+Ge.Range(1,gen))/2 + ...
                    randn(N,1)*(Ge.Range(2,gen)-Ge.Range(1,gen))/6;
            otherwise
                error(['Distribution type ',dist,' unknown']);
        end
        
        % Make sure the genome is within range
        Seq(:,gen) = max(min(Seq(:,gen),Ge.Range(2,gen)),Ge.Range(1,gen));
    end
    
    for s = 1:N
        % Check if the seq is valid
        [Res,Seq(s,:)] = Ge.CheckGenome(Seq(s,:));
        while Res{1} == 0
            % Create new random genome
            switch dist
                case {'uni','uniform'}
                    Seq(s,:) = Ge.Range(1,:)+rand(1,NGenes).*(Ge.Range(2,:)-Ge.Range(1,:));
                case {'norm','normal'}
                    Seq(s,:) = (Ge.Range(2,:)+Ge.Range(1,:))/2 + ...
                        randn(1,NGenes).*(Ge.Range(2,:)-Ge.Range(1,:))/6;
            end

            % Make sure the genome is within range
            Seq(s,:) = max(min(Seq(s,:),Ge.Range(2,:)),Ge.Range(1,:));
            [Res,Seq(s,:)] = Ge.CheckGenome(Seq(s,:));
        end
    end
end

