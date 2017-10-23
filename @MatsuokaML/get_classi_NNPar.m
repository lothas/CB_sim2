function seq = get_classi_NNPar(obj, NN, rand_seq)
%GET_CLASSI_NNPAR Use neural network to obtain new good seq
%   it recieve 'N' rand CPGs and classifi them to 'n-osc' and 'osc'
%   finally it chosses in random one good one.
    
    [genes_id] = obj.Gen.GetGenesId(obj.sample_genes);
    
    % transpose the array
    rand_seq = rand_seq';
    
    % prepare NN inputs:
    NN_in = rand_seq(genes_id,:);
    
    % check using the classifier which ones are oscillatory:
    check = NN(NN_in);
    [~,ind] = max(check);
    
    % find one "good" CPG:
    %   "good" CPG class is '2'
    rand_good_ind = randsample(find(ind==2),1);
    seq = rand_seq(:,rand_good_ind);
    
    % flip to get a raw seq vector:
    seq = seq';
end