function [seq_after_NN,NNoutputs] = train_and_test_NN(obj,caseNum,...
    architecture,method,testOn)
% this function test the NN for the wanted case.
% if we check the 'training' data, then we also plot the NN perf.

% *) 'caseNum' - the wanted case from the paper
% *) 'architecture' - row vector with the NN architecture
% *) 'method' - either: 'NN' - for neural networ
%                       'MoE' - for Mixture of Experts
% *) 'testOn' - wether to test on a random data, external data, or the same
% training data.
% 

%% NN training:
% get the names of the training parameters:
[Inputs_names_train,Targets_names_train] =...
    obj.check_NN_case(caseNum,'period');

% make NN inputs & outputs matrices:
[Inputs_train,Targets_train] = ...
    obj.prepare_NN_train_data(Inputs_names_train,Targets_names_train);

switch method
    case 'NN'
        % train NN:
        obj = obj.train_NN(architecture,Inputs_train,Targets_train);
    case 'MoE colaboration'
        numOfIteretions = 50;
        obj = obj.MoE_train_collaboration(numOfIteretions,...
            architecture,Inputs_train,Targets_train);
    case 'MoE hard'
        numOfIteretions = 50;
        obj = obj.MoE_train_Competetive('hardCompetetetive',numOfIteretions,...
            architecture,Inputs_train,Targets_train);
    case 'MoE soft'
        numOfIteretions = 50;
        obj = obj.MoE_train_Competetive('softCompetetetive',numOfIteretions,...
            architecture,Inputs_train,Targets_train);
    otherwise
        error('invalid method...');
end


%% NN testing:
[Inputs_names_test,Targets_names_test] = ...
    obj.check_NN_case(caseNum,'period_desired');

switch testOn
    case 'test_on_random_data'
        seq_test  = MML.Gen.RandSeq(100000);
    case 'test_on_training_data'
        seq_test  = vertcat(obj.results(obj.osc_ids).seq);
    case 'test_on_external_data'
        % choose 1000 random non-osc CPGs:
        id_rand = randsample();
        % TODO: finish to write this...
        
    otherwise
        error('invalid "testOn" case');
end

% acrivate the NN and change the right parameters according to the NN
% output:
[seq_after_NN,NNoutputs] = obj.prepare_NN_test_seq(seq_test,...
    Inputs_names_test,Targets_names_test,method);

%% plot:
net = obj.NN.net;
NNoutputs_on_train = net(Inputs_train);

targetsNum = size(Targets_names_train,2);

% plot histograms to compare the outputs to the inputs:
if strcmp(testOn,'test_on_training_data')
    
    if false % norm the outputs or not
        % norm the targets:
        temp = zeros(size(Targets_train));
        for i=1:targetsNum
            temp(i,:) = ...
                obj.norm_min_max(Targets_train(i,:),Targets_names_train{1,i});
        end
        Targets_train = temp;
        
        % norm the outputs::
        temp = zeros(size(NNoutputs_on_train));
        for i=1:targetsNum
            temp(i,:) = ...
                obj.norm_min_max(NNoutputs_on_train(i,:),...
                Targets_names_train{1,i});
        end
        NNoutputs_on_train = temp;
        
        % set the axis between '0' to '1':
        Axis = [0,1,0,1];
    else
        Axis = zeros(1,4);
        for i=1:targetsNum
            gen_id = strcmp(Targets_names_train{1,i},obj.seqOrder);
            Axis(1,2*i-1) = obj.MML.Gen.Range(1,gen_id)
            Axis(1,2*i) = obj.MML.Gen.Range(2,gen_id)
        end
    end
    
    switch targetsNum
        case 1
            figure;
            histogram(Targets_train,100,'Normalization','pdf'); hold on;
            histogram(NNoutputs_on_train,100,'Normalization','pdf');
            xlabel(Targets_names_train);
            legend('NNtargets','NNoutputs');
            grid minor;

        case 2
            
            xedges = linspace(Axis(1,1),Axis(1,2),100);
            yedges = linspace(Axis(1,3),Axis(1,4),100);
            
            figure;
            histogram2(Targets_train(1,:),Targets_train(2,:),...
                xedges,yedges,...
                'DisplayStyle','tile','ShowEmptyBins','on',...
                'Normalization','pdf');
            xlabel(Targets_names_train{1,1});
            ylabel(Targets_names_train{1,2});
            title('training data:');
            axis(Axis);
            
            figure;
            histogram2(NNoutputs_on_train(1,:),...
                NNoutputs_on_train(2,:),...
                xedges,yedges,...
                'DisplayStyle','tile','ShowEmptyBins','on',...
                'Normalization','pdf');
            xlabel(Targets_names_train{1,1});
            ylabel(Targets_names_train{1,2});
            title('NN output');
            axis(Axis);
    end

    % plot the outputs agains the inputs and calculate RMSE:
    for i=1:targetsNum
        
        % only plotting random 5000 samples:
        rand_ind = randsample(size(Targets_train,2),5000);
        
        RMSE = sqrt(immse(Targets_train(i,rand_ind),...
            NNoutputs_on_train(i,rand_ind)));
        disp(['the RMSE of ',Targets_names_train{1,i},' is: ',...
            num2str(RMSE)]);

        figure;
        scatter(Targets_train(i,rand_ind),...
            NNoutputs_on_train(i,rand_ind)); hold on;
        xlabel('targets');
        ylabel('outputs');
        title({'reggression graph: target vs. outputs',...
            Targets_names_train{1,i}});
        grid minor
        axis([0,1,0,1]);

    end
end
end
