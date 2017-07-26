
clc; clear all; close all;

% load CPGs results:
load('MatsRandomRes_4Neurons_4Paper_2.mat')

% get the CPGs periods
periods = horzcat(results(:).periods);

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68, 0.78];
MML.perLimOut = MML.perLim + [-0.08, 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% Filter CPG's where not both signals oscillating:
osc_ids = ~isnan(periods);
osc_ids = osc_ids(1,:) & osc_ids(2,:);

% Filter CPG's where the is a big difference between hip and ankle:
periods_ratios = (periods(1,:)./periods(2,:));
diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

% get Good CPG's
ids = osc_ids & diff_ids;

% get the mean of the hip and ankle periods:
periodsMean = mean(periods(:,ids),1);

% GET CPGs with period in range:
ids_des_period_before = ids & ((periods(1,:) > MML.perLimOut(1)) & (periods(1,:) < MML.perLimOut(2)));

results_osc = results(ids);     % oscillatory CPGs
results_rescaled = results_osc;

N = length(results_osc);

%% Running Rescaling:
tic
t_cur = tic;

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
% for i=1:N
    disp(['at sim #',num2str(i)]);
    
    % check if the period is not in range:
    if ( (periodsMean(1,i) < MML.perLimOut(1)) || (periodsMean(1,i) > MML.perLimOut(2)) )
        seq = results_osc(i).seq;
        
        % Select new random period within desired range
        des_period = MML.perLim(1) + rand()*(MML.perLim(2)-MML.perLim(1));
        
        % Scale Tr, Ta to obtain desired period
        ratio = des_period/periodsMean(1,i);
        seq(1) = seq(1)*ratio; % change 'tau' in the sequence
        
        if seq(1) < MML.Gen.Range(1,1) || seq(1) > MML.Gen.Range(2,1)
            warning('Genetic sequence out of bounds, using bounded tau gene')
            % Bound tau gene
            seq(1) = min(max(seq(1), MML.Gen.Range(1,1)), MML.Gen.Range(2,1));
        end
        
        [out, sim, signal] = MML.runSim(seq);
        
            % Prepare output:
        % Parameters
        results_rescaled(i).seq = seq;
        results_rescaled(i).b = sim.Con.beta;
        results_rescaled(i).c = sim.Con.Amp0;
        results_rescaled(i).Worig = sim.Con.wex;
        results_rescaled(i).W = sim.Con.W;
        results_rescaled(i).Tr = sim.Con.tau;
        results_rescaled(i).Ta = sim.Con.tav;
        results_rescaled(i).x0 = out.x0;

        % Results- caculate perdiods using different methods:
        results_rescaled(i).periods = out.periods;

        results_rescaled(i).pos_work = out.pos_work;
        results_rescaled(i).neg_work = out.neg_work;
        results_rescaled(i).perError1 = out.perError1;
        results_rescaled(i).perOK1 = out.perOK1;
        results_rescaled(i).perError2 = out.perError2;
        results_rescaled(i).perOK2 = out.perOK2;
        results_rescaled(i).neuronActive = out.neuronActive;
        results_rescaled(i).neuronOsc = out.neuronOsc;
    end

end 
disp('sim end...');

t_elapsed = toc(t_cur);
avg_sim_time = t_elapsed/N;
disp(['avg sim time is ',num2str(avg_sim_time),' [sec]']);

periods2 = horzcat(results_rescaled(:).periods);
% Filter CPG's where not both signals oscillating:
osc_ids2 = ~isnan(periods2);
osc_ids2 = osc_ids2(1,:) & osc_ids2(2,:);

% Filter CPG's where the is a big difference between hip and ankle:
periods_ratios2 = (periods2(1,:)./periods2(2,:));
diff_ids2 = (periods_ratios2 >  0.85) & (periods_ratios2 <  1.15); 

ids2 = osc_ids2 & diff_ids2;

ids_des_period_after = ids2 & ...
     ((periods2(1,:) > MML.perLimOut(1)) & (periods2(1,:) < MML.perLimOut(2)));

disp(['the number of CPG in range before rescaling is: ',...
    num2str(sum(ids_des_period_before))]);
disp(['the number of CPG in range after rescaling is: ',...
    num2str(sum(ids_des_period_after))]);


save('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims_2.mat',...
    'results_rescaled','results_osc');


%% Rescale a random CPG:

i = randsample(sum(ids),1);

disp(['at CPG #',num2str(i)]);

% check if the period is not in range:
if ( (periodsMean(1,i) < MML.perLimOut(1)) || (periodsMean(1,i) > MML.perLimOut(2)) )
    disp(['The CPG need rescaling! the period is: ',...
        num2str(periodsMean(1,i))]);
    
    seq = results_osc(i).seq;

    % Select new random period within desired range
    des_period = MML.perLim(1) + rand()*(MML.perLim(2)-MML.perLim(1));

    % Scale Tr, Ta to obtain desired period
    ratio = des_period/periodsMean(1,i);
    seq(1) = seq(1)*ratio; % change 'tau' in the sequence

    if seq(1) < MML.Gen.Range(1,1) || seq(1) > MML.Gen.Range(2,1)
        disp(['    tau is ',seq(1)]);
        warning('Genetic sequence out of bounds, using bounded tau gene')
        % Bound tau gene
        seq(1) = min(max(seq(1), MML.Gen.Range(1,1)), MML.Gen.Range(2,1));
    end

    [out, sim, signal] = MML.runSim(seq);

    % Results- caculate perdiods using different methods:
    periods_cur = out.periods;
    
    disp(['    The period of the CPG after is ',num2str(periods_cur(1)),...
        ' and ',num2str(periods_cur(2))]);
else
    disp(['The period of the CPG is ',num2str(periodsMean(1,i)),...
        ' the CPG is in range']);

end
