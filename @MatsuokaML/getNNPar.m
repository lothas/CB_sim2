function seq = getNNPar(obj, NN, seq, period)
%GETNNPAR Use neural network to obtain new values of tau
    in = obj.getNNin(seq, period);
    seq(obj.Gen.GetGenesId(obj.target_genes)) = NN(in)';
%     NNout = NN(in);
%     
%     Tr = NNout(1);
end