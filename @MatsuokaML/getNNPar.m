function [Tr, Ta] = getNNPar(obj, NN, c, W, period)
%GETNNPAR Use neural network to obtain new values of tau
    in = [c;
         W(~logical(eye(size(W))));
         period];

    % Normalize X
    in = bsxfun(@rdivide, bsxfun(@minus, in, ...
        obj.normParams(:,1)), obj.normParams(:,2));
    NNout = NN(in);
    Tr = NNout(1);
    Ta = obj.TrTaRatio*NNout(1);
end