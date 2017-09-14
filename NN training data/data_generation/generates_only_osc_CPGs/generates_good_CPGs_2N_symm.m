
close all; clc; clear all;

%%
genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 0;   % Number of hip torques
maxAnkle = 10;   % Max ankle torque
maxHip = 10;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
Mw = 10*ones(1,(2*N-1)*2*N);
mw = 0*Mw;
% %     % 2neuron symmetric specific range%%

% % % % Large b Large W
% Mw = 10*ones(1,(2*N-1)*2*N);
% mw = 0*Mw;
% Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                         1,        1 ,          2,            0 };
% Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
%            0.25  ,    10,      [10,10],                        10,    0.001 , [0.2,0.2]}; % Max
% 

%        % % % Large b Narrow W
% Mw = 10*ones(1,(2*N-1)*2*N);
% mw = 0*Mw;
% Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                         1,        1 ,          2,            0 };
% Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
%            0.25  ,    10,      [10,10],                         5,    0.001 , [0.2,0.2]}; % Max

       % % % Narrow b Narrow W
Mw = 5*ones(1,(2*N-1)*2*N);
mw = 0*Mw;
Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,            2,                         1,        1 ,          2,            0 };
Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
           0.25  ,     5,      [10,10],                         5,    0.001 , [0.2,0.2]}; % Max

       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%%
% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2;

%% Train data: (get good ones)
N = 1000; % the number of samples
good_CPGs_num = 0;
results = [];

wanted_num_CPGs = 100000;

disp('start with the sim:');

while good_CPGs_num < wanted_num_CPGs
    rand_seq = MML.Gen.RandSeq(N);
%     results_temp(N) = [];
    parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
    % for i=1:N
%         disp(['at sim #',num2str(i)]);
        [out, ~, signal] = MML.runSim(rand_seq(i,:));
            % Prepare output:
        % Parameters
        results_temp(i).seq = rand_seq(i,:);

        % Results- caculate perdiods using different methods:
        results_temp(i).periods = out.periods;  

        results_temp(i).pos_work = out.pos_work;
        results_temp(i).neg_work = out.neg_work;
        results_temp(i).perError1 = out.perError1;
        results_temp(i).perOK1 = out.perOK1;
        results_temp(i).perError2 = out.perError2;
        results_temp(i).perOK2 = out.perOK2;
        results_temp(i).neuronActive = out.neuronActive;
        results_temp(i).neuronOsc = out.neuronOsc;
    end 

    % get the periods
    periods = horzcat(results_temp(:).periods);
    % get oscillatory CPGs ids:
    osc_ids = ~isnan(periods);

    % get neural oscillation check:
    neuronOsc = (vertcat(results_temp(:).neuronOsc))';
    osc_check_ids = true(1,N);
    for k=1:size(neuronOsc,1) % run on all neurons
        % mark as "good" if we have at least one osc neuron
        osc_check_ids = osc_check_ids & ~neuronOsc(k,:);
    end

    % get perOK1 check:
    perOK1 = (vertcat(results_temp(:).perOK1))';

    good_ids = osc_ids & ~osc_check_ids & perOK1;

    results = [results,results_temp(good_ids)];
    good_CPGs_num = length(results);

    disp(['so far we have ',num2str(good_CPGs_num),' CPGs'])

    clear results_temp periods rand_seq osc_check_ids neuronOsc
    clear osc_check_ids good_ids osc_ids

end

disp('sim end...');


% header = sprintf('tau ratio is equal to 12 \n');
% header = [header,sprintf('data is for 2N symmetric case \n')];
% header = [header,sprintf('seq Order: \n')];
% header = [header,sprintf('"tau","b","c","NR","a" \n')];
% header = [header,sprintf('"b" in range (0.2,2.5) \n')];
% header = [header,sprintf('"a" in range (0,5) \n')];
% 
% 
% save('MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_only_osc_3.mat',...
%     'results','header','MML');

%% Train data: (get uniformly dist of 'b')
N = 1000; % the number of samples
good_CPGs_num = 0;
results = [];
results_temp = [];
results_store = [];

wanted_num_CPGs = 100;
b_ranges = linspace(0.2001,2.4999,100);

disp('start with the sim:');

for j=2:length(b_ranges)
    good_CPGs_num = 0;
    disp(['at bin #',num2str(j-1)]);
    while good_CPGs_num < wanted_num_CPGs
        MML.Gen.Range(1,2) = b_ranges(1,j-1);
        MML.Gen.Range(2,2) = b_ranges(1,j);
        rand_seq = MML.Gen.RandSeq(N);
    %     results_temp(N) = [];
        parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
        % for i=1:N
    %         disp(['at sim #',num2str(i)]);
            [out, ~, signal] = MML.runSim(rand_seq(i,:));
                % Prepare output:
            % Parameters
            results_temp(i).seq = rand_seq(i,:);

            % Results- caculate perdiods using different methods:
            results_temp(i).periods = out.periods;  

            results_temp(i).pos_work = out.pos_work;
            results_temp(i).neg_work = out.neg_work;
            results_temp(i).perError1 = out.perError1;
            results_temp(i).perOK1 = out.perOK1;
            results_temp(i).perError2 = out.perError2;
            results_temp(i).perOK2 = out.perOK2;
            results_temp(i).neuronActive = out.neuronActive;
            results_temp(i).neuronOsc = out.neuronOsc;
        end 

        % get the periods
        periods = horzcat(results_temp(:).periods);
        % get oscillatory CPGs ids:
        osc_ids = ~isnan(periods);

        % get neural oscillation check:
        neuronOsc = (vertcat(results_temp(:).neuronOsc))';
        osc_check_ids = true(1,N);
        for k=1:size(neuronOsc,1) % run on all neurons
            % mark as "good" if we have at least one osc neuron
            osc_check_ids = osc_check_ids & ~neuronOsc(k,:);
        end

        % get perOK1 check:
        perOK1 = (vertcat(results_temp(:).perOK1))';

        good_ids = osc_ids & ~osc_check_ids & perOK1;

        results_store = [results_store,results_temp(good_ids)];
        good_CPGs_num = length(results_store);

        clear results_temp periods rand_seq osc_check_ids neuronOsc
        clear osc_check_ids good_ids osc_ids

    end
    results = [results,...
        results_store(randsample(1:length(results_store),wanted_num_CPGs))];
    disp(['so far we have ',num2str(length(results)),' CPGs']);
end
disp('sim end...');


% header = sprintf('tau ratio is equal to 12 \n');
% header = [header,sprintf('data is for 2N symmetric case \n')];
% header = [header,sprintf('seq Order: \n')];
% header = [header,sprintf('"tau","b","c","NR","a" \n')];
% header = [header,sprintf('"b" in range (0.2,2.5) \n')];
% header = [header,sprintf('"a" in range (0,5) \n')];
% 
% 
% save('MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_only_osc_3.mat',...
%     'results','header','MML');
%% plot CPG output:
close all;

n = randsample(1:length(results),1);

[out, ~, signal] = MML.runSim(results(n).seq);
figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
title({'X_i over time',...
    ['id #',num2str(n),...
    '    periods: ',...
    num2str(out.periods)]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
clear signal
