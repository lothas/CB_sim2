function [a, b, c, Worig, W, Tr, Ta] = getRandPar(obj)
%GETRANDPAR Generate new random parameters for the CPG
    a = 2;
    b = 4;
    c = rand(obj.nNeurons,1)*10;

    W = rand(obj.nNeurons, obj.nNeurons);
%     W = W*2+(1-W)*-1;
    W = 3*W - 1;
    W(logical(eye(size(W)))) = 0;

    blks = cell(obj.nNeurons/2, 1);
    blks(:) = {~eye(2)};
    Worig = W + a*blkdiag(blks{:});

    % Normalize W by c
    W = ((Worig*diag(1./c))'*diag(c))';
    % Keep weights in range of +/- 20
    W = min(max(W, -20), 20); %#ok<UDIM>

    Tr = 0.01+rand()*0.36; % 0.01*exp(5.5*rand());
    Ta = obj.TrTaRatio*Tr; % 0.01*exp(5.5*rand());
end


