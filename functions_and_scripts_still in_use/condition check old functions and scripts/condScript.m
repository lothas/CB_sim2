
clc; close all; clear all;

MML = MatsuokaML();
MML.tStep = 0.01;
MML.tEnd = 15;
%% Load varibles:
load('MatsRandomRes_4_11_2016.mat','nSims','results','periods')
load('neuronActive_4_11_2016.mat');

%% building the condition matrix (or loading an existing one):

filename = 'condMat_4_11_2016.mat';
if ~exist(filename,'file')
    condMat = checkCondMat(nSims,results,periods,neuronActive);
    save(filename,'condMat');
else
    load(filename);
end


%% define the flags:
% condMat = [ A_positive ,cond0, cond1, cond2, cond2UC, cond2_Sim, cond2UC_Sim, has_period];
A_positive = condMat(:,1);
cond0 = condMat(:,2);
neuron_i_active = [condMat(:,3),condMat(:,4),condMat(:,5),condMat(:,6)];
% at least one active neuron present
al_1_active_N_fromCond = (condMat(:,3) | condMat(:,4) | condMat(:,5) | condMat(:,6));
al_1_active_N_fromSim = (neuronActive(:,1) | neuronActive(:,2) | neuronActive(:,3) | neuronActive(:,4));
% by checking the necessary condition for active neuron:
oscSol_cond2 = condMat(:,7);
oscSol_cond2UC = condMat(:,8);
% by checking the simulation results for active neuron:
oscSol_cond2_Sim = condMat(:,9);
oscSol_cond2UC_Sim = condMat(:,10);
% seq has a detected period:
hasPeriod = condMat(:,end);
% all N are active:
all_N_are_active_fromCond = (condMat(:,3) & condMat(:,4) & condMat(:,5) & condMat(:,6));
all_N_are_active_fromSim = (neuronActive(:,1) & neuronActive(:,2) & neuronActive(:,3) & neuronActive(:,4));
% How many active neurons:
howManyActiveN_fromCond = sum([condMat(:,3),condMat(:,4),condMat(:,5),condMat(:,6)],2);
howManyActiveN_fromSim = sum([neuronActive(:,1),neuronActive(:,2),neuronActive(:,3),neuronActive(:,4)],2);
% 3 N are active:
three_N_are_active_fromCond = (howManyActiveN_fromCond == 3);
three_N_are_active_fromSim = (howManyActiveN_fromSim == 3);
% 2 N are active:
two_N_are_active_fromCond = (howManyActiveN_fromCond == 2);
two_N_are_active_fromSim = (howManyActiveN_fromSim == 2);
% 1 N are active:
one_N_are_active_fromCond = (howManyActiveN_fromCond == 2);
one_N_are_active_fromSim = (howManyActiveN_fromSim == 2);

%% For cond2 results in the paper
clc
disp('truth table: for the "real" cond2: ');
tT_1st_cond = A_positive & oscSol_cond2 ;
tT_2nd_cond = A_positive & condMat(:,end);
rowsNames = {['cond2',' true:'],['cond2',' false:']};
columnNames = {['Oscillating',' true:'],['Oscillating',' false:']};
tT_11 = sum(tT_1st_cond & tT_2nd_cond);
tT_12 = sum(tT_1st_cond & ~tT_2nd_cond);
tT_21 = sum(~tT_1st_cond & tT_2nd_cond);
tT_22 = sum(~tT_1st_cond & ~tT_2nd_cond);
truthTable = [0             ,columnNames(1,1),columnNames(1,2);
              rowsNames(1,1),tT_11           ,tT_12 ;
              rowsNames(1,2),tT_21           ,tT_22];
disp(truthTable)
clear tT_11 tT_12 tT_21 tT_22 rowsNames columnNames tT_1st_cond tT_2nd_cond 

disp('Same table with percetage');
tT_1st_cond = A_positive & oscSol_cond2 ;
tT_2nd_cond = A_positive & condMat(:,end);
rowsNames = {['cond2',' true:'],['cond2',' false:']};
columnNames = {['Oscillating',' true:'],['Oscillating',' false:']};
tT_11 = sum(tT_1st_cond & tT_2nd_cond)*100/sum(tT_1st_cond);
tT_12 = sum(tT_1st_cond & ~tT_2nd_cond)*100/sum(tT_1st_cond);
tT_21 = sum(~tT_1st_cond & tT_2nd_cond)*100/sum(~tT_1st_cond);
tT_22 = sum(~tT_1st_cond & ~tT_2nd_cond)*100/sum(~tT_1st_cond);
truthTable = [0             ,columnNames(1,1),columnNames(1,2);
              rowsNames(1,1),tT_11           ,tT_12 ;
              rowsNames(1,2),tT_21           ,tT_22];
disp(truthTable)
disp(['cond2 true: ',num2str(sum(tT_1st_cond)*100/sum(A_positive))]);
clear tT_11 tT_12 tT_21 tT_22 rowsNames columnNames tT_1st_cond tT_2nd_cond 


%% For cond1 results in the paper
clc
disp('truth table: for cond1: ');
neuron1cond =  A_positive & condMat(:,3) ;
neuron2cond =  A_positive & condMat(:,4) ;
neuron3cond =  A_positive & condMat(:,5) ;
neuron4cond =  A_positive & condMat(:,6) ;
neuron1act =  A_positive & neuronActive(:,1);
neuron2act =  A_positive & neuronActive(:,2);
neuron3act =  A_positive & neuronActive(:,3);
neuron4act =  A_positive & neuronActive(:,4);
rowsNames = {['cond1',' true:'],['cond1',' false:']};
columnNames = {['N_i active',' true:'],['N_i active',' false:']};
tT_11 = sum(neuron1cond & neuron1act) + sum(neuron2cond & neuron2act) + ...
    sum(neuron3cond & neuron3act) + sum(neuron4cond & neuron4act);
tT_12 = sum(neuron1cond & ~neuron1act) + sum(neuron2cond & ~neuron2act) + ...
    sum(neuron3cond & ~neuron3act) + sum(neuron4cond & ~neuron4act);
tT_21 = sum(~neuron1cond & neuron1act) + sum(~neuron2cond & neuron2act) + ...
    sum(~neuron3cond & neuron3act) + sum(~neuron4cond & neuron4act);
tT_22 = sum(~neuron1cond & ~neuron1act) + sum(~neuron2cond & ~neuron2act) + ...
    sum(~neuron3cond & ~neuron3act) + sum(~neuron4cond & ~neuron4act);
truthTable = [0             ,columnNames(1,1),columnNames(1,2);
              rowsNames(1,1),tT_11           ,tT_12 ;
              rowsNames(1,2),tT_21           ,tT_22];
disp(truthTable)

disp('Same table with percetage');
tT_11a = tT_11*100/(sum(neuron1cond)+sum(neuron2cond)+sum(neuron3cond)+sum(neuron4cond));
tT_12a = tT_12*100/(sum(neuron1cond)+sum(neuron2cond)+sum(neuron3cond)+sum(neuron4cond));
tT_21a = tT_21*100/(sum(~neuron1cond)+sum(~neuron2cond)+sum(~neuron3cond)+sum(~neuron4cond));
tT_22a = tT_22*100/(sum(~neuron1cond)+sum(~neuron2cond)+sum(~neuron3cond)+sum(~neuron4cond));
truthTable = [0             ,columnNames(1,1),columnNames(1,2);
              rowsNames(1,1),tT_11a          ,tT_12a ;
              rowsNames(1,2),tT_21a          ,tT_22a];
disp(truthTable)
clear tT_11a tT_12a tT_21a tT_22a rowsNames columnNames
clear tT_11 tT_12 tT_21 tT_22
clear neuron1cond neuron2cond neuron3cond neuron4cond
clear neuron1act neuron2act neuron3act neuron4act


%% checking all of the examples "without Period" to see whether it's true:           
figure;
N = find(tT_1st_cond & ~tT_2nd_cond);
for j=1:78
    sequence = results(N(j)).seq;

    [~, ~, signal] = runSim(MML, sequence);
    
    if mod(j,20)==0
        figure;
    end
    subplot(4,5,mod(j,20)+1)
    plot(signal.T,signal.y); grid minor;
    title(['#seq = ',num2str(N(j)), ' and period = ', num2str(periods(1,N(j)))]);
    xlabel('time[sec]'); ylabel('y_i');

end

clear tT_11 tT_12 tT_21 tT_22 rowsNames columnNames tT_1st_cond tT_2nd_cond 