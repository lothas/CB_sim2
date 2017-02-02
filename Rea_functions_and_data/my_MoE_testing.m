function [R_squar,errMSE,netOut,belongToExpert] = my_MoE_testing(NNinputs,NNtargets,expertsNN,...
    gateNet,GraphicsFlag,competetiveFlag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: as from 1/2/2017 this function is obsolete!
% DO NOT USE. use instead "my_MoE_testNet"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function train a "Mixture of Experts" newural networks

% inputs:
% 1) NNinputs - the inputs for the neural networks,
%               matrix of (#inputs x #samples) dimention
% 2) NNtargets - the tarhets for the neural networks,
%                matrix of (#targets x #samples) dimention
% 3) expertsNN - cell array containing the experts and their performance.
% 4) gateNet - the Gate NN.
% 5) competetiveFlag - if 'competetive'- than the clustering is being done
%                       by highest P takes all. if not competetive than
%                       each sample as a chance to go to each cluster based
%                       on the probability from the gate network

% outputs:
% 1) R_squar - R^2 after rechecking the training data
% 2) errMSE - MSE after rechecking the training data
% 3) netOut - output from MoE
% 4) belongToExpert - which expert's cluster the sample is belong.

expertCount = size(expertsNN,2);

% check the final clustering after training:
belongToExpert = zeros(1,size(NNinputs,2));
gateOut = gateNet(NNinputs);
if competetiveFlag
    [~,belongToExpert] = max(gateOut,[],1); % indecies of points to #clusters
else
   % cluster by probability and not by highest chance
   for k=1:size(gateOut,2)
        acumulativeProb = tril(ones(expertCount,expertCount))*gateOut(:,k);
        belongToExpert(1,k) = find(acumulativeProb > rand(1),1);
   end
end

if size(NNinputs,2) > 1 % many samples, performance analisys
    cluster_i_train_ind = cell(1,expertCount);
    targ_train = [];
    netOut = [];
    if  GraphicsFlag 
        colors = rand(expertCount,3);
        legendNames = cell(1,expertCount);
        for j=1:expertCount
            legendNames{1,j} = ['#',num2str(j),' expert'];
        end
        figure; hold on
    end

    for j=1:expertCount
        cluster_i_train_ind{1,j} = find(gateOut_inx == j);
        tempNet = expertsNN{1,j};
        out_train_temp = tempNet(NNinputs(:,cluster_i_train_ind{1,j}));
        targ_temp = NNtargets(:,cluster_i_train_ind{1,j});
        targ_train = [targ_train,targ_temp];
        netOut = [netOut,out_train_temp];

        if GraphicsFlag
            h = plot(targ_temp,out_train_temp,'Color',colors(j,:),'LineStyle','none');
            h.Marker = 'o';
        end
        clear targ_temp out_train_temp
    end

    if GraphicsFlag
        hold off;
        xlabel('targets'); ylabel('ouput'); legend(legendNames);
    end

    % show performance:
    [errMSE,R_squar] = NN_perf_calc(targ_train,netOut,1,GraphicsFlag );
    
else % single sample, return MoE output
    tempNet = expertsNN{1,belongToExpert};
    netOut = tempNet(NNinputs(:,1));
    errMSE = [];
    R_squar = [];
end

end

