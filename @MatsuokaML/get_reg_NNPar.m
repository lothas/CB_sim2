function seq = get_reg_NNPar(obj, NN, seq, period)
%GET_REG_NNPAR Use neural network to obtain new values of 'tau' or 'b'
    in = obj.getNNin(seq, period);
    seq(obj.Gen.GetGenesId(obj.target_genes)) = NN(in)';

end