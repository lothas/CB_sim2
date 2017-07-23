function in = getNNin(obj, seq, period)
%GETNNPAR Use neural network to obtain new values of tau
    % Get tonic inputs and connection weights from genetic sequence
    if obj.norm_weights
        genes1 = obj.Gen.GetGenes(seq, obj.sample_genes1); % weights
        amp = obj.Gen.GetGenes(seq, {'amp'}); % weights
        for r = 1:obj.nNeurons
            range = (obj.nNeurons-1)*(r-1)+1:(obj.nNeurons-1)*r;
            this_amp = amp(r);
            amps = amp; amps(r) = [];
            genes1(range) = genes1(range)*this_amp./amps;
        end
%         genes1 = min(max(genes1, -10),100); % Bound the weights
        genes2 = obj.Gen.GetGenes(seq, obj.sample_genes2);
        genes = [genes1, genes2];
    else
        genes = obj.Gen.GetGenes(seq, obj.sample_genes);
    end
    
%     in = [genes';  1/period];
    in = [genes';  period];

    % Normalize X
    if ~isempty(obj.normParams)
        in = bsxfun(@rdivide, bsxfun(@minus, in, ...
            obj.normParams(:,1)), obj.normParams(:,2));
    end
end