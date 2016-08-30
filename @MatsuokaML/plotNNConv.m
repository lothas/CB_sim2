function varargout = plotNNConv(obj, rndData, sclData, NNData, type)
%PLOTNNCONV Summary of this function goes here
%   Detailed explanation goes here

    % Build data:
    if type == 1
        % Architecture labels
        xval = 1:1+NNData.nArch;
        xname = cell(1,1+NNData.nArch);
        xname{1} = 'rand';
        for i = 1:NNData.nArch
            xname{1+i} = num2str(NNData.architectures{i});
        end
    else
        % Number of data points labels
        sampScale = 10000.0;
        zero = 2*NNData.sampleSizes(1) - NNData.sampleSizes(2);
        xval = [zero, NNData.sampleSizes];
        xname = cell(1,1+NNData.nSampleSizes);
        xname{1} = 'rand';
        for i = 1:NNData.nSampleSizes
            xname{1+i} = num2str(NNData.sampleSizes(i)/sampScale);
        end
    end
    nPoints = length(xval);

    % Percentage of simulations that converged in general and within the
    % desired period range
    NNSamples = size(NNData.sampPerf,2);
    genConv = zeros(1,nPoints);
    perConv = zeros(1,nPoints);
    genConvSc = zeros(1,nPoints);
    perConvSc = zeros(1,nPoints);
    genConv(1) = length(rndData.id_conv)/rndData.nSims;
    genConvSc(1) = length(sclData.id_conv)/rndData.nSims;
    perConv(1) = length(rndData.id_per)/rndData.nSims;
    perConvSc(1) = length(sclData.id_per)/rndData.nSims;
    for i = 1:nPoints-1
        % Find scaled results that converged
        ids = find(~isnan(NNData.sampPerf(i,:)));
        genConv(1+i) = length(ids)/NNSamples;
        % Find how many converged within the desired period
        pers = NNData.sampPerf(i,ids);
        within = pers >= obj.perLimOut(1) & pers <= obj.perLimOut(2);
        perConv(1+i) = sum(within)/NNSamples;

        % Find scaled results that converged
        scIds = find(~isnan(NNData.sampPerfSc(i,:)));
        genConvSc(1+i) = length(scIds)/NNSamples;
        % Find how many converged within the desired period
        scPers = NNData.sampPerfSc(i,scIds);
        within = scPers >= obj.perLimOut(1) & scPers <= obj.perLimOut(2);
        perConvSc(1+i) = sum(within)/NNSamples;
    end

    if nargout == 1
        % Prepare output data
        dataOut.xval = xval;
        dataOut.xname = xname;
        dataOut.genConv = genConv;
        dataOut.genConvSc = genConvSc;
        dataOut.perConv = perConv;
        dataOut.perConvSc = perConvSc;
        varargout = {dataOut};
    else
        % Plot results
        figure
        plot(xval,genConv,'*-')
        hold on
        plot(xval,genConvSc,'s-')
        plot(xval,perConv,'o-')
        plot(xval,perConvSc,'^-')
        legend('Sim converged', 'Sim converged (with scaling)', ...
            'Sim conv. in range', 'Sim conv. in range (with scaling)', ...
            'Location','Southeast')
        set(gca,'XTick',xval,'XTickLabel',xname)
        ylabel('Percentage convergence')
        if type == 1
            xlabel('Architectures')
            title('Comparison of approaches')
        else
            xlabel(['Number of training samples [x', ...
                int2str(sampScale),'] (log)'])
            title(['Performance of [', ...
                num2str(NNData.architecture),'] net'])
            set(gca,'xscale','log')
        end

        figure
        % plot(xval(2:end),genConv(2:end),'*-')
        hold on
        plot(xval(2:end),genConvSc(2:end),'s-')
        % plot(xval(2:end),perConv(2:end),'o-')
        plot(xval(2:end),perConvSc(2:end),'^-')
        legend('Sim converged', 'Sim conv. in range', ...
            'Location','Southeast')
        set(gca,'XTick',xval(2:end),'XTickLabel',xname(2:end))
        title('NN+Scaling performance comparison')
        ylabel('Percentage convergence')
        if type == 1
            xlabel('NN architecture')
        else
            xlabel(['Number of training samples [x', ...
                int2str(sampScale),'] (log)'])
            set(gca,'xscale','log')
        end
        
        varargout = {};
    end
end

