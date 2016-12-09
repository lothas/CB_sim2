function pSeq = permSeq( obj, seq, new_order )
%PERMSEQ Summary of this function goes here
%   Detailed explanation goes here

doPermC = 0;
if nargin<3
    % Allowed permutations for 2 hip and 2 ankle neurons
    new_order = [1,2,3,4;
                 2,1,3,4;
                 1,2,4,3;
                 2,1,4,3];
             
    % If we only permutate W
    if ~doPermC
        new_order = perms(1:4);
    end
end

% Rebuild W and C
Wids = obj.Gen.GetGenesId({'weights'});
Cids = obj.Gen.GetGenesId({'amp'});
nNeurons = length(Cids);

C = seq(Cids)';
W = zeros(nNeurons, nNeurons);
k = 1;
for i = 1:nNeurons
    for j = 1:nNeurons
        if i == j
            continue
        end
        W(i,j) = seq(Wids(k));
        k = k+1;
    end
end

% Permutate
nPerms = size(new_order,1);
pSeq = repmat(seq, nPerms, 1);
for i = 1:nPerms
    permC = C(new_order(i,:));
    permW = W(new_order(i,:),new_order(i,:));

    % Re-encode sequence
    if doPermC
        pSeq(i,Cids) = permC';
    end
    permW = reshape(permW',1,[]); permW(permW == 0) = [];
    pSeq(i,Wids) = permW;
end

end

