function [Inputs_names,Targets_names] = check_NN_case(obj,caseNum,perORfreq)
% This function check which of the 9 different NN cases we want to load
% then, make 'inputs_names' and 'target_names' accordinly.


% Inputs: *)'caseNum': |#num  |NN input            |   NN output    |
%                      | '1'  |W_ij_hat            |tau_r           |
%                      | '2'  |ci, W_ij_hat        |tau_r           |
%                      | '3'  |b, W_ij_hat         |tau_r           |
%                      | '4'  |b, ci, W_ij_hat     |tau_r           |
%                      | '5'  |W_ij_hat            |b               |
%                      | '6'  |ci, W_ij_hat        |b               |
%                      | '7'  |tau_r, W_ij_hat     |b               |
%                      | '8'  |tau_r, ci, W_ij_hat |b               |
%                      | '9'  |W_ij_hat            |tau_r, b        |
%                      | '10' |ci, W_ij_hat        |tau_r, b        |
%                      | '11' |tau_r               |b               |
%          *) 'perORfreq' - whter to choose to use the period or the
%                           frequency

% outputs: *) 'parametersCells' - NN inputs names
%          *) 'targetCells' - NN outputs names

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% NOTE: Jonathan is arranging his NN differently
% FIRST: the CPG parameters
% THAN after: the Period
% e.g. [seq,period]
% % % % % % % % % % % % % % % % % % % % % % !!!!!!!!!!!


switch perORfreq
    case 'period'
        theta_s0 = {'period'};
    case 'freq'
        theta_s0 = {'freq'};
    case 'period_desired'
        theta_s0 = {'period_desired'};
    case 'freq_desired'
        theta_s0 = {'freq_desired'};
    otherwise
        error('you can only choose between the period or the frequency!');
end

W_ij = {'w_{12}','w_{13}','w_{14}',...
        'w_{21}','w_{23}','w_{24}',...
        'w_{31}','w_{32}','w_{34}',...
        'w_{41}','w_{42}','w_{43}'};
ci = {'c_1','c_2','c_3','c_4'};
    
switch caseNum
    case 1  % NN input outputs: (option 1)
        tempCell = {[theta_s0,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'tau'};
    case 2  % NN input outputs: (option 2)
        tempCell = {[theta_s0,ci,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'tau'};
    case 3  % NN input outputs: (option 3)
        tempCell = {[theta_s0,{'b'},W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'tau'};
    case 4  % NN input outputs: (option 4)
        tempCell = {[theta_s0,{'b'},ci,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'tau'};
    case 5  % NN input outputs: (option 5)
        tempCell = {[theta_s0,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'b'};
    case 6  % NN input outputs: (option 6)
        tempCell = {[theta_s0,ci,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'b'};
    case 7  % NN input outputs: (option 7)
        tempCell = {[theta_s0,{'tau'},W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'b'};
    case 8  % NN input outputs: (option 8)
        tempCell = {[theta_s0,{'tau'},ci,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'b'};
    case 9  % NN input outputs: (option 9)
        tempCell = {[theta_s0,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'tau','b'};
    case 10  % NN input outputs: (option 10)
        tempCell = {[theta_s0,ci,W_ij]};
        parametersCells = tempCell{1,1};
        targetCells = {'tau','b'};
    case 11
        tempCell = {[theta_s0,{'tau'}]};
        parametersCells = tempCell{1,1};
        targetCells = {'b'};
        
    otherwise
        error('wrong case number');
end

Inputs_names = parametersCells;
Targets_names = targetCells;

end

