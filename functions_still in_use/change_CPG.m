function [results_new] = ...
    change_CPG(MML,caseNum,seqOrder,theta_S1_new,results_old)
% this function make a new CPG from the theta_S1 and theta_S2(from NN)

% Inputs: *)'caseNum': |#num  |NN input (theta_S2) | NN output (theta_S1)|
%                      | '1'  |W_ij_hat            |tau_r                |
%                      | '2'  |ci, W_ij_hat        |tau_r                |
%                      | '3'  |b, W_ij_hat         |tau_r                |
%                      | '4'  |b, ci, W_ij_hat     |tau_r                |
%                      | '5'  |W_ij_hat            |b                    |
%                      | '6'  |ci, W_ij_hat        |b                    |
%                      | '7'  |tau_r, W_ij_hat     |b                    |
%                      | '8'  |tau_r, ci, W_ij_hat |b                    |
%                      | '9'  |W_ij_hat            |tau_r, b             |
%                      | '10' |ci, W_ij_hat        |tau_r, b             |
%          *) 'seqOrder' - the order of the parameters in the gene
%          *) 'theta_S1_new' - the new parameters from the NN
%          *) 'results_old' - the old CPG parameters
%          *) 'MML' - MML class for the simulation

% outputs: *) 'results_new' - new CPG to run
%         

seq = vertcat(results_old(:).seq);

switch caseNum
    case {1,2,3,4}
        change_index = strcmp('tau',seqOrder);
        seq(:,change_index) = theta_S1_new'; % theta_S1_new = tau_new
    case {5,6,7,8}
        change_index = strcmp('b',seqOrder);
        seq(:,change_index) = theta_S1_new'; % theta_S1_new = b
    case {9,10}
        change_index = strcmp('tau',seqOrder) | strcmp('b',seqOrder);
        theta_S1_new=theta_S1_new';
        seq(:,change_index) = theta_S1_new; % theta_S1_new =tau
end

N = size(seq,1);
if 1 % 'if' for code folding purposes
    [out, ~, ~] = MML.runSim(seq(N,:));
        % Prepare output:
    % Parameters
    results_new(N).seq = seq(N,:);
    results_new(N).x0 = out.x0;
    % Results
    results_new(N).periods = out.periods;
    results_new(N).pos_work = out.pos_work;
    results_new(N).neg_work = out.neg_work;
    results_new(N).perError1 = out.perError1;
    results_new(N).perOK1 = out.perOK1;
    results_new(N).perError2 = out.perError2;
    results_new(N).perOK2 = out.perOK2;
    results_new(N).neuronActive = out.neuronActive;
    results_new(N).neuronOsc = out.neuronOsc;
end    

parfor i=1:(N-1)
    [out, ~, ~] = MML.runSim(seq(i,:));
    
        % Prepare output:
    % Parameters
    results_new(i).seq = seq(i,:);
    results_new(i).x0 = out.x0;

    % Results
    results_new(i).periods = out.periods;
    results_new(i).pos_work = out.pos_work;
    results_new(i).neg_work = out.neg_work;
    results_new(i).perError1 = out.perError1;
    results_new(i).perOK1 = out.perOK1;
    results_new(i).perError2 = out.perError2;
    results_new(i).perOK2 = out.perOK2;
    results_new(i).neuronActive = out.neuronActive;
    results_new(i).neuronOsc = out.neuronOsc;

end
results_new = results_new';

end

