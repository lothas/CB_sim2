function [err, r, g] = paper_MoE_test(obj,sampl, targ, ExpertsWeights, gateWeights,GraphicsFlag)
  % based on code from the internet: https://goker.wordpress.com/2011/07/01/mixture-of-experts/
  
    sampl = sampl';
    targ = targ';
    outputCount = size(targ, 2);
    expertCount = size(gateWeights,1);
    sampleCount = size(sampl, 1);
    
    % add bias
    sampl = [ones(sampleCount, 1), sampl];
    
    ExpertOut = zeros(outputCount, expertCount, sampleCount);
    
    for j=1:outputCount
        ExpertOut(j, :, :) = (ExpertsWeights(:,:,j)*sampl');
    end
    
    g = (exp(gateWeights*sampl'));        
    g = g ./ repmat(sum(g, 1), expertCount, 1);
    y = zeros(sampleCount, outputCount);
    
    for i=1:sampleCount
        y(i, :) = ( ExpertOut(:,:,i) * g(:,i) );        
    end

    err = sum((y - targ).^2) ./ sampleCount;
    r = y';
    
    if GraphicsFlag
        figure;
        subplot(2,1,1);
        h = plot(targ,y,'LineStyle','none'); grid minor;
        h.Marker = 'o';
        xlabel('targets'); ylabel('outputs');
        title('regression graph: Targets over NNoutputs');
        
        legendNames = cell(1,expertCount);
        for j=1:expertCount
            legendNames{1,j} = ['#',num2str(j),' expert'];
        end
    
        subplot(2,1,2);
        bar(g','stacked'); xlabel('#sample');
        legend(legendNames);
        
    end
 
end