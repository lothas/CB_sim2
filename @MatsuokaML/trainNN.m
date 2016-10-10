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
%     getRandFuncHandle = @obj.Gen.RandSeq;
    genomeObj = obj.Gen;
    getNNFuncHandle = @obj.getNNPar;
    runSimFuncHandle = @obj.runSim;
    perLimMin = obj.perLim(1);
    perLimMax = obj.perLim(2);
    parfor j = 1:NNSamples
        % Setup tau_r, tau_a, c, W and feedback gains using genome
%         seq = getRandFuncHandle(); % Get random genetic sequence
        seq = genomeObj.RandSeq(); %#ok<PFBNS> % Get random genetic sequence

%         if ~any(strcmp(genomeObj.Keys(1,:),'beta'))
%             % Set random b
%             if rand()>0.7
%                 beta = min(max(0.6+0.1*randn(),0.2),0.8);
%             else
%                 beta = min(max(2.5+randn(),0.8),8);
%             end
%         end
        
        % Get Tr estimates from NN
        seqNN = getNNFuncHandle(net, seq, desPeriod(j));
        [res, seqNN] = genomeObj.CheckGenome(seqNN);
        if (res{1})
            seq = seqNN;
        else
            warning(['Genetic sequence #', int2str(j), ...
                ' out of bounds, keeping original sequence'])
        end
        
        % Run simulation
        [out, ~, ~] = runSimFuncHandle(seq);
%         [out, ~, ~] = runSimFuncHandle(seq, beta);
        
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
            seq(1) = seq(1)*ratio;
            if seq(1) < genomeObj.Range(1,1) || ...
                    seq(1) > genomeObj.Range(2,1)
                warning(['Genetic sequence #', int2str(j), ...
                    ' out of bounds, using bounded tau gene'])
                % Bound tau gene
                seq(1) = min(max(seq(1), genomeObj.Range(1,1)), ...
                             genomeObj.Range(2,1));
            end

            % Run simulation
            [out, ~, ~] = runSimFuncHandle(seq);
%             [out, ~, ~] = runSimFuncHandle(seq, beta);
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

