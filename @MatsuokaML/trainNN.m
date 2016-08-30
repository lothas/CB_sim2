function [net, tr, netPerf, desPeriod, sampPerf, sampPerfSc] = ...
    trainNN(obj, samples, targets, architecture, NNSamples)
%TRAINNN Train a NN and calculate its performance
    netPerf = zeros(1, 4);    % Array to store NN performance
    sampPerf = zeros(1, NNSamples);
    sampPerfSc = zeros(1, NNSamples);

    % Create and train the NN
    net = feedforwardnet(architecture);
    [net, tr] = train(net, samples, targets);
        
    % Calculate mean squared estimation error
    estTargets = net(samples);
    netPerf(1) = sum(sum((estTargets - targets).^2))/size(targets,2);
                
    % Generate array of desired period outputs for the simulations
    desPeriod = obj.perLim(1) + ...
                 rand(1, NNSamples)*(obj.perLim(2)-obj.perLim(1));
             
    % Run adaptive simulations using the NN output and see how many
    % converge to the right period
    getRandFuncHandle = @obj.getRandPar;
    getNNFuncHandle = @obj.getNNPar;
    runSimFuncHandle = @obj.sim;
    perLimMin = obj.perLim(1);
    perLimMax = obj.perLim(2);
    parfor j = 1:NNSamples
        % Setup oscillators' parameters
        [~, b, c, ~, W, ~, ~] = getRandFuncHandle();
        % Get Tr, Ta estimates from NN
        [Tr, Ta] = getNNFuncHandle(net, c, W, desPeriod(j));
        % Run simulation
        [out, ~] = runSimFuncHandle(b, c, W, Tr, Ta);
        
        % Save resulting period
        if any(isnan(out.periods))
            sampPerf(j) = NaN;
        else
            sampPerf(j) = max(out.periods);
        end
        
        % Did period converge to the desired range?
        if sampPerf(j) < perLimMin || sampPerf(j) > perLimMax
            % Try rescaling time
            ratio = desPeriod(j)/sampPerf(j);
            Tr = Tr*ratio;
            Ta = 5*Tr;
            [out, ~] = runSimFuncHandle(b, c, W, Tr, Ta);
            % Save resulting period
            if any(isnan(out.periods))
                sampPerfSc(j) = NaN;
            else
                sampPerfSc(j) = max(out.periods);
            end
        else
            sampPerfSc(j) = sampPerf(j);
        end
    end
    
    % How many converged?
    id_conv = find(~isnan(sampPerf));
    netPerf(2) = length(id_conv)/NNSamples;
    % How many converged to the desired range?
    id_per = find(sampPerf >= obj.perLimOut(1) ...
        & sampPerf <= obj.perLimOut(2));
    netPerf(3) = length(id_per)/NNSamples;
    % How close was the period to the desired period?
    sampDiff = sampPerf - desPeriod;
    netPerf(4) = sum(sampDiff(id_conv).^2)/NNSamples;
end

