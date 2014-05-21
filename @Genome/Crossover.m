function [ BabySeqs, Res ] = Crossover( Ge, Seq1, Seq2, Method )
%CROSSOVER Creates 2 new genome sequences by crossing 2 given sequences
%   Crossover uses a 2-point crossover to mix the given sequences
%   and create two new ones, which are checked before being returned
    if nargin<4
        Method = '2point';
    end
    
    N = length(Seq1);
    Bits = zeros(2,N);
    
    switch Method
        case {'1point','one-point','one point'}
            % Choose 1 crossover point at random
            Point = randsample(1:N-1,1);
            Bits(1,:) = (1:N)<=Point;  % Child 1
            Bits(2,:) = ~Bits(1,:);    % Child 2
        case {'2point','two-point','two point'}
            % Choose 2 crossover points at random
            Point = sort(randsample(1:N-1,2));
            Bits(1,:) = (1:N)<=Point(1) | (1:N)>Point(2);   % Child 1
            Bits(2,:) = ~Bits(1,:);                         % Child 2
        case {'uni','uniform'}
            % Crossover with 50% probability
            bit = 1;
            for b = 1:N
                Bits(1,b) = bit;       % Child 1
                if rand(1)>0.5
                    bit = ~bit;
                end
            end
            Bits(2,:) = ~Bits(1,:);    % Child 2
        otherwise
            error(['Crossover method ',Method,' unknown']);
    end
    
    BabySeqs = repmat(Seq1,2,1).*Bits+repmat(Seq2,2,1).*~Bits;
    Res = cell(2,2);
    [Res(1,:), BabySeqs(1,:)] = Ge.CheckGenome(BabySeqs(1,:));
    [Res(2,:), BabySeqs(2,:)] = Ge.CheckGenome(BabySeqs(2,:));
end

