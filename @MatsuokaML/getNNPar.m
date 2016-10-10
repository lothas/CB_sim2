function seq = getNNPar(obj, NN, seq, period)
%GETNNPAR Use neural network to obtain new values of tau
    % Get tonic inputs and connection weights from genetic sequence
    genes = obj.Gen.GetGenes(seq, obj.sample_genes);
    in = [genes';  period];

    % Normalize X
    in = bsxfun(@rdivide, bsxfun(@minus, in, ...
        obj.normParams(:,1)), obj.normParams(:,2));
    seq(obj.Gen.GetGenesId(obj.target_genes)) = NN(in)';
%     NNout = NN(in);
%     
%     Tr = NNout(1);
end