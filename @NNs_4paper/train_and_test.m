function obj = train_and_test(obj,Inputs_names,Targets_names,...
    architecture,method,graphFlag)
% this function test the NN for the wanted case.
% if we check the 'training' data, then we also plot the NN perf.

% *) 'architecture' - row vector with the NN architecture
% *) 'method' - either: 'NN' - for neural networ
%                       'MoE' - for Mixture of Experts
% *) 'testOn' - wether to test on a random data, external data, or the same
% training data.
% *) 'graphFlag' - plot the graphs or not

%% NN training:

% make NN inputs & outputs matrices:
[Inputs_train,Targets_train] = ...
    obj.prepare_NN_train_data(Inputs_names,Targets_names);

switch method
    case 'NN'
        % train NN:
        obj.NN = NeuralNetwork(architecture,Inputs_train,Targets_train);
        obj.NN = obj.NN.train_NN();
        targ_train = Targets_train(:,obj.NN.net_train_perf.trainInd);
%         obj.NN.net_train_function = 'trainrp';
    case 'MoE colaboration'
        obj.MoE = MoE(Inputs_train,Targets_train,'collaboration');
        obj.MoE.numOfIteretions = 5;
        obj.MoE = obj.MoE.MoE_train();
        targ_train = Targets_train(:,obj.MoE.train_Ind);
    case 'MoE hard'
        obj.MoE = MoE(Inputs_train,Targets_train,'hardCompetetetive');
        obj.MoE.numOfIteretions = 5;
        obj.MoE = obj.MoE.MoE_train();
        targ_train = Targets_train(:,obj.MoE.train_Ind);
    case 'MoE soft'
        obj.MoE = MoE(Inputs_train,Targets_train,'softCompetetetive');
        obj.MoE.numOfIteretions = 5;
        obj.MoE = obj.MoE.MoE_train();
        targ_train = Targets_train(:,obj.MoE.train_Ind);

    otherwise
        error('invalid method...');
end

%% plot: test NN on the training data

targetsNum = size(Targets_names,2);

switch method
    case 'NN'
        in = Inputs_train(:,obj.NN.net_train_perf.trainInd);
        % TODO: change to give only the train group...
    case {'MoE colaboration','MoE hard','MoE soft'}
        in = Inputs_train(:,obj.MoE.train_Ind);
end
NNoutputs_on_train = obj.apply_net(in,method);

% plot histograms to compare the outputs to the inputs:
if graphFlag
    
    if false % norm the outputs or not
        % norm the targets:
        temp = zeros(size(targ_train));
        for i=1:targetsNum
            temp(i,:) = ...
                obj.norm_min_max(targ_train(i,:),Targets_names{1,i});
        end
        Targets_train = temp;
        
        % norm the outputs::
        temp = zeros(size(NNoutputs_on_train));
        for i=1:targetsNum
            temp(i,:) = ...
                obj.norm_min_max(NNoutputs_on_train(i,:),...
                Targets_names{1,i});
        end
        NNoutputs_on_train = temp;
        
        % set the axis between '0' to '1':
        Axis = [0,1,0,1];
    else
        Axis = zeros(1,4);
        for i=1:targetsNum
            gen_id = strcmp(Targets_names{1,i},obj.seqOrder);
            Axis(1,2*i-1) = obj.MML.Gen.Range(1,gen_id);
            Axis(1,2*i) = obj.MML.Gen.Range(2,gen_id);
        end
    end
    
    switch targetsNum
        case 1
            figure;
            histogram(targ_train,100,'Normalization','pdf'); hold on;
            histogram(NNoutputs_on_train,100,'Normalization','pdf');
            xlabel(Targets_names);
            legend('NNtargets','NNoutputs');
            grid minor;

        case 2
            
            xedges = linspace(Axis(1,1),Axis(1,2),100);
            yedges = linspace(Axis(1,3),Axis(1,4),100);
            
            figure;
            histogram2(targ_train(1,:),targ_train(2,:),...
                xedges,yedges,...
                'DisplayStyle','tile','ShowEmptyBins','on',...
                'Normalization','pdf');
            xlabel(Targets_names{1,1});
            ylabel(Targets_names{1,2});
            title('training data:');
            axis(Axis);
            
            figure;
            histogram2(NNoutputs_on_train(1,:),...
                NNoutputs_on_train(2,:),...
                xedges,yedges,...
                'DisplayStyle','tile','ShowEmptyBins','on',...
                'Normalization','pdf');
            xlabel(Targets_names{1,1});
            ylabel(Targets_names{1,2});
            title('NN output');
            axis(Axis);
    end

    % plot the outputs agains the inputs and calculate RMSE:
    for i=1:targetsNum
        
        % only plotting random 5000 samples:
        rand_ind = randsample(size(targ_train,2),5000);

        figure;
        scatter(targ_train(i,rand_ind),...
            NNoutputs_on_train(i,rand_ind)); hold on;
        xlabel('targets');
        ylabel('outputs');
        title({'reggression graph: target vs. outputs',...
            Targets_names{1,i}});
        grid minor
        axis([0,1,0,1]);

    end
end
end
